"""FastAPI entrypoint: validate, place, and deliver jplace results."""

from __future__ import annotations

import json
import uuid

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.responses import Response
from pydantic import BaseModel, Field, conint

from .community import (build_sample_mass, geometry_of,
                        pairwise_distances)
from .jplacefmt import build_jplace
from .likelihood import DirectedMessages, obs_vector
from .placement import has_base_evidence, place_query
from .tree import from_phylotree, orient_and_serialize, renumber_edges
from .validation import (MAX_QUERIES, MAX_REFS, MIN_QUERIES, MIN_REFS,
                         SubmissionError, cross_validate, parse_tree,
                         validate_alignment)

app = FastAPI(title="eDNA placement backend", version="1.1.0")

JOBS: dict[str, dict] = {}

MIN_SAMPLES, MAX_SAMPLES = 2, 10


class PlacementRequest(BaseModel):
    reference_fasta: str
    query_fasta: str
    newick: str


class CompareSample(BaseModel):
    sample_id: str = Field(min_length=1)
    job_id: str = Field(min_length=1)
    counts: dict[str, conint(ge=0)]


class CompareRequest(BaseModel):
    samples: list[CompareSample] = Field(min_length=MIN_SAMPLES,
                                         max_length=MAX_SAMPLES)


def run_placement(req: PlacementRequest) -> dict:
    problems: list[str] = []
    ref = validate_alignment(req.reference_fasta, "reference",
                             MIN_REFS, MAX_REFS, problems)
    qry = validate_alignment(req.query_fasta, "query",
                             MIN_QUERIES, MAX_QUERIES, problems)
    tree = parse_tree(req.newick, problems)
    if ref is not None and qry is not None and tree is not None:
        cross_validate(ref, qry, tree, problems)
    if problems:
        raise SubmissionError(problems)

    g = from_phylotree(tree)
    renumber_edges(g)
    tree_str = orient_and_serialize(g)

    name_to_node = {name: n for n, name in g.leaf_name.items()}
    leaf_obs = {}
    for sid, states in zip(ref.ids, ref.states):
        leaf_obs[name_to_node[sid]] = np.array([obs_vector(s) for s in states])
    dm = DirectedMessages(g, leaf_obs)

    placements_by_query: dict[str, list] = {}
    unplaceable: dict[str, str] = {}
    for sid, states in zip(qry.ids, qry.states):
        if not has_base_evidence(states):
            unplaceable[sid] = ("no valid base evidence: every column is a gap "
                                "or fully ambiguous, so the sequence carries no "
                                "phylogenetic signal and cannot be placed")
            continue
        query_obs = np.array([obs_vector(s) for s in states])
        placements_by_query[sid] = place_query(dm, g, query_obs)

    jplace = build_jplace(tree_str, placements_by_query)
    return {"jplace": jplace, "unplaceable": unplaceable,
            "edge_count": len(g.edges)}


@app.post("/place", status_code=200)
def place(req: PlacementRequest):
    try:
        result = run_placement(req)
    except SubmissionError as exc:
        raise HTTPException(status_code=422,
                            detail={"rejected": True, "problems": exc.problems})
    job_id = uuid.uuid4().hex[:12]
    JOBS[job_id] = result
    summary = {
        "job_id": job_id,
        "jplace_url": "/place/" + job_id + "/jplace",
        "edge_count": result["edge_count"],
        "unplaceable": result["unplaceable"],
        "placements": {
            p["n"][0]: [
                dict(zip(result["jplace"]["fields"], row)) for row in p["p"]
            ]
            for p in result["jplace"]["placements"]
        },
    }
    return summary


@app.get("/place/{job_id}/jplace")
def download_jplace(job_id: str):
    if job_id not in JOBS:
        raise HTTPException(status_code=404, detail="unknown job_id")
    payload = json.dumps(JOBS[job_id]["jplace"], indent=2)
    return Response(
        content=payload,
        media_type="application/json",
        headers={"Content-Disposition":
                 "attachment; filename=placements_" + job_id + ".jplace"},
    )


def run_compare(req: CompareRequest) -> dict:
    problems: list[str] = []
    samples = req.samples
    seen_ids: set[str] = set()
    for pos, smp in enumerate(samples):
        if smp.sample_id in seen_ids:
            problems.append("sample at position " + str(pos) + ": duplicate "
                            "sample_id '" + smp.sample_id + "'")
        seen_ids.add(smp.sample_id)
        if not smp.counts:
            problems.append("sample '" + smp.sample_id + "': counts map is empty")
        job = JOBS.get(smp.job_id)
        if job is None:
            problems.append("sample '" + smp.sample_id + "': unknown job_id '" +
                            smp.job_id + "'")
            continue
        known = set(job["unplaceable"]) | {
            p["n"][0] for p in job["jplace"]["placements"]}
        for qid in smp.counts:
            if qid not in known:
                problems.append("sample '" + smp.sample_id + "': unknown query "
                                "id '" + qid + "' for job '" + smp.job_id + "'")

    known_jobs = [JOBS[smp.job_id] for smp in samples if smp.job_id in JOBS]
    if known_jobs:
        ref_tree = known_jobs[0]["jplace"]["tree"]
        for smp in samples[1:]:
            if smp.job_id not in JOBS:
                continue
            job = JOBS[smp.job_id]
            if job["jplace"]["tree"] != ref_tree:
                problems.append("sample '" + smp.sample_id + "': reference tree "
                                "(annotated Newick) is not identical to the "
                                "first sample's job; community comparison only "
                                "supports tasks placed on the same tree")
    if problems:
        raise SubmissionError(problems)

    built = [build_sample_mass(smp.sample_id, smp.counts, JOBS[smp.job_id])
             for smp in samples]
    for smp in built:
        if smp.effective_total <= 0.0:
            problems.append("sample '" + smp.sample_id + "': effective total "
                            "is zero (no readable, placeable counts)")
    if problems:
        raise SubmissionError(problems)

    geom = geometry_of(JOBS[samples[0].job_id])
    matrix, pairs = pairwise_distances(built, geom)
    return {
        "samples": [smp.sample_id for smp in samples],
        "distance_matrix": matrix,
        "pairs": pairs,
        "excluded": [
            {"sample_id": smp.sample_id, "items": smp.excluded}
            for smp in built if smp.excluded
        ],
        "effective_totals": {
            smp.sample_id: smp.effective_total for smp in built},
    }


@app.post("/compare", status_code=200)
def compare(req: CompareRequest):
    try:
        return run_compare(req)
    except SubmissionError as exc:
        raise HTTPException(status_code=422,
                            detail={"rejected": True,
                                    "problems": exc.problems})
