"""Self-test: run a placement in-process and sanity-check the results."""

import json
import pathlib

from fastapi.testclient import TestClient

from app.main import app

ROOT = pathlib.Path(__file__).resolve().parent.parent
EX = ROOT / "examples"


def main() -> None:
    client = TestClient(app)
    payload = {
        "reference_fasta": (EX / "reference.fasta").read_text(),
        "query_fasta": (EX / "queries.fasta").read_text(),
        "newick": (EX / "tree.nwk").read_text(),
    }
    resp = client.post("/place", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert "q3" in data["unplaceable"], data
    for qid in ("q1", "q2"):
        places = data["placements"][qid]
        assert len(places) == data["edge_count"] == 7
        weights = [p["like_weight_ratio"] for p in places]
        assert abs(sum(weights) - 1.0) < 1e-9, weights
        lls = [p["likelihood"] for p in places]
        assert lls == sorted(lls, reverse=True)
        for p in places:
            assert 0.0 <= p["pendant_length"] <= 2.0
            assert p["distal_length"] >= 0.0
    # q1 matches refA exactly: its best edge should touch refA's pendant edge
    best_q1 = data["placements"]["q1"][0]
    print("q1 best edge:", best_q1)

    dl = client.get("/place/" + data["job_id"] + "/jplace")
    assert dl.status_code == 200
    doc = json.loads(dl.content)
    assert doc["version"] == 3
    assert doc["fields"] == ["edge_num", "likelihood", "like_weight_ratio",
                             "distal_length", "pendant_length"]
    assert "{" in doc["tree"] and "}" in doc["tree"]
    assert len(doc["placements"]) == 2

    # rejection path: illegal character must be located and reject the batch
    bad = dict(payload)
    bad["query_fasta"] = ">qX\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGX\n"
    resp = client.post("/place", json=bad)
    assert resp.status_code == 422
    problems = resp.json()["detail"]["problems"]
    assert any("qX" in p and "column 40" in p for p in problems), problems
    print("rejection OK:", problems[0])
    # DNA must not accept RNA base U
    rna = dict(payload)
    rna["query_fasta"] = ">qU\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGU\n"
    resp = client.post("/place", json=rna)
    assert resp.status_code == 422
    up = resp.json()["detail"]["problems"]
    assert any("qU" in p and "column 40" in p and "'U'" in p for p in up), up
    print("U rejection OK")

    # leaf ids containing commas must survive Newick export quoted
    comma_tree = ("((refA:0.1,refB:0.12):0.08,(refC:0.11,refD:0.09):0.07,"
                  "'sp,E':0.2);")
    comma_ref = (EX / "reference.fasta").read_text().replace(">refE\n",
                                                             ">sp,E\n")
    cresp = client.post("/place", json={
        "reference_fasta": comma_ref,
        "query_fasta": (EX / "queries.fasta").read_text(),
        "newick": comma_tree})
    assert cresp.status_code == 200, cresp.text
    doc = client.get("/place/" + cresp.json()["job_id"] + "/jplace").json()
    assert "'sp,E'" in doc["tree"], doc["tree"]
    from app.jplacetree import parse_jplace_tree
    gback = parse_jplace_tree(doc["tree"])
    assert len(gback.leaf_name) == 5
    assert "sp,E" in gback.leaf_name.values()
    print("comma leaf quoting OK")

    # ---- community comparison ----
    comp = {"samples": [
        {"sample_id": "site1", "job_id": data["job_id"],
         "counts": {"q1": 5, "q2": 3, "q3": 2}},
        {"sample_id": "site2", "job_id": data["job_id"],
         "counts": {"q1": 1, "q2": 8}},
        {"sample_id": "site3", "job_id": data["job_id"],
         "counts": {"q1": 4}},
    ]}
    resp = client.post("/compare", json=comp)
    assert resp.status_code == 200, resp.text
    out = resp.json()
    assert out["samples"] == ["site1", "site2", "site3"]
    mat = out["distance_matrix"]
    assert len(mat) == 3 and all(len(row) == 3 for row in mat)
    assert all(mat[i][i] == 0.0 for i in range(3))
    assert all(abs(mat[i][j] - mat[j][i]) < 1e-12 for i in range(3)
               for j in range(3))
    assert out["effective_totals"] == {"site1": 8.0, "site2": 9.0,
                                       "site3": 4.0}
    excl = {e["sample_id"]: e["items"] for e in out["excluded"]}
    assert excl["site1"][0]["id"] == "q3"
    assert excl["site1"][0]["count"] == 2
    assert "cannot be placed" in excl["site1"][0]["reason"]
    assert len(out["pairs"]) == 3
    for p in out["pairs"]:
        contribs = p["edge_contributions"]
        assert len(contribs) == data["edge_count"]
        assert [c["edge_num"] for c in contribs] == list(
            range(data["edge_count"]))
        assert all(c["contribution"] >= 0.0 for c in contribs)
        assert abs(sum(c["contribution"] for c in contribs) -
                   p["distance"]) < 1e-9
    assert abs(mat[0][1] - out["pairs"][0]["distance"]) < 1e-12
    print("compare OK; matrix:",
          [[round(v, 4) for v in row] for row in mat])

    def expect_422(label, body):
        rr = client.post("/compare", json=body)
        assert rr.status_code == 422, (label, rr.status_code, rr.text)
        detail = rr.json()["detail"]
        return detail["problems"] if isinstance(detail, dict) else [str(detail)]

    one = {"sample_id": "a", "job_id": data["job_id"], "counts": {"q1": 1}}
    two_good = {"sample_id": "b", "job_id": data["job_id"],
                "counts": {"q1": 1}}
    rr = client.post("/compare", json={"samples": [one]})
    assert rr.status_code == 422, rr.text
    pr = expect_422("unknown job", {"samples": [
        {"sample_id": "a", "job_id": "deadbeef", "counts": {"q1": 1}},
        two_good]})
    assert any("unknown job_id" in p for p in pr), pr
    pr = expect_422("unknown query", {"samples": [
        {"sample_id": "a", "job_id": data["job_id"], "counts": {"nope": 1}},
        two_good]})
    assert any("unknown query" in p and "nope" in p for p in pr), pr
    expect_422("negative count", {"samples": [
        {"sample_id": "a", "job_id": data["job_id"], "counts": {"q1": -2}},
        two_good]})
    pr = expect_422("duplicate sample", {"samples": [
        {"sample_id": "a", "job_id": data["job_id"], "counts": {"q1": 1}},
        {"sample_id": "a", "job_id": data["job_id"], "counts": {"q2": 1}}]})
    assert any("duplicate" in p for p in pr)
    pr = expect_422("zero effective total", {"samples": [
        {"sample_id": "a", "job_id": data["job_id"], "counts": {"q3": 3}},
        two_good]})
    assert any("effective total" in p for p in pr), pr
    pr = expect_422("tree mismatch", {"samples": [
        {"sample_id": "a", "job_id": data["job_id"], "counts": {"q1": 1}},
        {"sample_id": "b", "job_id": cresp.json()["job_id"],
         "counts": {"q1": 1}}]})
    assert any("reference tree" in p for p in pr), pr
    print("compare rejections OK")
    print("SELFTEST PASSED")


if __name__ == "__main__":
    main()
