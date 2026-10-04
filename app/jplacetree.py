"""Parse the jplace 'tree' string back into oriented tree geometry.

A jplace tree is Newick with an integer edge-number annotation in braces
after every branch length, e.g. "(A:0.1{3},B:0.2{4}):0.05{2};".  The
parser keeps the parent->child orientation written in the string, so the
recovered Edge.u/Edge.v match distal_length semantics in the placement
records (distal = towards the child node v).
"""

from __future__ import annotations

from .tree import Edge, TreeGraph


class JplaceTreeError(ValueError):
    """Malformed jplace tree string."""


def parse_jplace_tree(text: str) -> TreeGraph:
    return _Parser(text.strip()).parse()


class _Parser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.pos = 0
        self.g = TreeGraph()
        self._raw_edges: list[tuple[int, int, float, int]] = []

    def parse(self) -> TreeGraph:
        self._clade(None)
        self._skip_ws()
        if self.pos >= len(self.text) or self.text[self.pos] != ";":
            raise JplaceTreeError("expected terminating ';'")
        self.pos += 1
        self._skip_ws()
        if self.pos != len(self.text):
            raise JplaceTreeError("trailing characters after ';'")
        nums = sorted(num for _, _, _, num in self._raw_edges)
        if nums != list(range(len(self._raw_edges))):
            raise JplaceTreeError(
                "edge numbers must be a permutation of 0.." +
                str(len(self._raw_edges) - 1))
        ordered = sorted(self._raw_edges, key=lambda t: t[3])
        self.g.edges = [Edge(u, v, length) for u, v, length, _ in ordered]
        self.g.adj = {n: [] for n in self.g.adj}
        for eid, e in enumerate(self.g.edges):
            self.g.adj[e.u].append((e.v, eid))
            self.g.adj[e.v].append((e.u, eid))
        return self.g

    def _skip_ws(self) -> None:
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1

    def _peek(self) -> str:
        self._skip_ws()
        return self.text[self.pos] if self.pos < len(self.text) else ""

    def _clade(self, parent: int | None) -> int:
        if self._peek() != "(":
            raise JplaceTreeError("expected '(' at position " + str(self.pos))
        self.pos += 1
        node = self.g.add_node()
        children: list[int] = []
        while True:
            ch = self._peek()
            if ch == "(":
                children.append(self._clade(node))
            elif ch and ch not in ",);":
                children.append(self._leaf(node))
            else:
                raise JplaceTreeError("unexpected character '" + ch + "'")
            sep = self._peek()
            if sep == ",":
                self.pos += 1
                continue
            if sep == ")":
                self.pos += 1
                break
            raise JplaceTreeError("expected ',' or ')' at position " +
                                  str(self.pos))
        self._branch(node, parent)
        return node

    def _leaf(self, parent: int) -> int:
        node = self.g.add_node()
        self.g.leaf_name[node] = self._label()
        self._branch(node, parent)
        return node

    def _label(self) -> str:
        self._skip_ws()
        if self.pos < len(self.text) and self.text[self.pos] == "'":
            self.pos += 1
            chars: list[str] = []
            while True:
                if self.pos >= len(self.text):
                    raise JplaceTreeError("unterminated quoted label")
                ch = self.text[self.pos]
                self.pos += 1
                if ch == "'":
                    if self.pos < len(self.text) and self.text[self.pos] == "'":
                        chars.append("'")
                        self.pos += 1
                        continue
                    break
                chars.append(ch)
            return "".join(chars)
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in ":,);":
            self.pos += 1
        name = self.text[start:self.pos].strip()
        if not name:
            raise JplaceTreeError("unnamed leaf at position " + str(start))
        return name

    def _branch(self, node: int, parent: int | None) -> None:
        """Parse ':length{edge_num}' for the branch above node."""
        self._skip_ws()
        if parent is None:
            if self.pos < len(self.text) and self.text[self.pos] == ":":
                raise JplaceTreeError("root branch length not allowed")
            return
        if self._peek() != ":":
            raise JplaceTreeError("missing branch length at position " +
                                  str(self.pos))
        self.pos += 1
        start = self.pos
        while self.pos < len(self.text) and (self.text[self.pos].isdigit()
                                            or self.text[self.pos] in ".eE+-"):
            self.pos += 1
        token = self.text[start:self.pos]
        try:
            length = float(token)
        except ValueError as exc:
            raise JplaceTreeError("bad branch length '" + token + "'") from exc
        self._skip_ws()
        if self._peek() != "{":
            raise JplaceTreeError("missing edge number annotation at position " +
                                  str(self.pos))
        self.pos += 1
        nstart = self.pos
        while self.pos < len(self.text) and self.text[self.pos].isdigit():
            self.pos += 1
        num_token = self.text[nstart:self.pos]
        if not num_token:
            raise JplaceTreeError("empty edge number at position " + str(nstart))
        if self._peek() != "}":
            raise JplaceTreeError("unclosed edge number at position " +
                                  str(self.pos))
        self.pos += 1
        self._raw_edges.append((parent, node, length, int(num_token)))
