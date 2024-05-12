#include <iostream>
#include <vector>
#include <unordered_set>
#include "NumCpp.hpp"
#include <limits>
using namespace std;

string tt(const nc::NdArray<double>& t) {
    return to_string(t[0]) + " " + to_string(t[1]);
}


class DelaunayTriangulation {
private:
    unordered_map<int, unordered_map<int, int>> edges;
    unordered_map<int, unordered_set<int>> edge_to_faces;
    unordered_map<int, vector<int>> face_to_vertex;
    nc::NdArray<double> vertices;

    int edgeId = 0;
    int faceId = 0;
public:
    DelaunayTriangulation(const nc::NdArray<double>& points) : vertices(points) {
        int idx_p0 = highestRight();
        {
            nc::NdArray<double> p = leftAndRightMost(idx_p0);
            vertices.row(idx_p0).print();
            p.print();
            vertices = nc::vstack({vertices, p});
        }
        int idx_p__1 = vertices.shape().rows-2;
        int idx_p__2 = vertices.shape().rows-1;

        this->addTriangle(idx_p0, idx_p__1, idx_p__2);

        for (int i = 0; i < points.shape().rows; i++) {
            if (i == idx_p0) continue;
            auto res = this->containInTriangle(i);
            if (res.second.first == -1 and res.second.second == 0) {
                // Inside
                auto vertexes = face_to_vertex[res.first];
                eraseTriangle(res.first);

                int f1 = addTriangle(i, vertexes[0], vertexes[1]);
                int f2 = addTriangle(i, vertexes[1], vertexes[2]);
                int f3 = addTriangle(i, vertexes[2], vertexes[0]);

                legalizeEdge(i, {vertexes[0], vertexes[1]}, f1);
                legalizeEdge(i, {vertexes[1], vertexes[2]}, f2);
                legalizeEdge(i, {vertexes[2], vertexes[0]}, f3);
            }
            else {
                // on edge
                int i_ = res.second.first;
                int j_ = res.second.second;
                auto faces = edge_to_faces[edges[i_][j_]];
                int k = complement(face_to_vertex[*faces.begin()], i_, j_);
                eraseTriangle(*faces.begin());
                int fk_i = addTriangle(k, i, i_);
                int fk_j = addTriangle(k, i, j_);
                int l = -1, fl_i, fl_j;
                if (faces.size() > 1) {
                    l = complement(face_to_vertex[*next(faces.begin())], i_, j_);
                    eraseTriangle(*next(faces.begin()));
                    fl_i = addTriangle(i, l, i_);
                    fl_j = addTriangle(i, l, j_);
                }

                if (l != -1) {
                    legalizeEdge(i, {i_, l}, fl_i);
                    legalizeEdge(i, {l, j_}, fl_j);
                }
                legalizeEdge(i, {j_, k}, fk_j);
                legalizeEdge(i, {k, i_}, fk_i);
            }
        }

        vector<int> face_to_erase;
        for (auto [k, face]: face_to_vertex) {
            auto it = std::find_if(face.begin(), face.end(), [&](const auto &item) {
                return item == idx_p__2 or item == idx_p__1;
            });
            if (it != face.end()) {
                face_to_erase.emplace_back(k);
            }
        }

        for (auto e: face_to_erase) {
            eraseTriangle(e);
        }
    }

    int highestRight() {
        auto val = -nc::constants::inf;
        auto idx = -1;
        for (int i = 0; i < vertices.shape().rows; i++){
            if (vertices(i, 1) > val) {
                val = vertices(i, 1);
                idx = i;
            }
            else if (vertices(i, 1) == val and vertices(i, 0) > vertices(idx, 0)) {
                idx = i;
            }
        }
        return idx;
    }

    nc::NdArray<double> leftAndRightMost(int p) {
        nc::NdArray<double> axis = {1, 0};
        auto mn_val = nc::constants::inf;
        int mn_idx = -1;
        auto mx_val = -nc::constants::inf;
        int mx_idx = -1;
        for (int i = 0; i < vertices.shape().rows; i++) {
            if (i == p) continue;
            auto v = vertices.row(i) - vertices.row(p);
            auto tmp = (nc::dot(v, axis) / nc::norm(v))[0];
            if (tmp < mn_val) {
                mn_val = tmp;
                mn_idx = i;
            }
            if (tmp > mx_val) {
                mx_val = tmp;
                mx_idx = i;
            }
        }
        nc::NdArray<double> left = vertices.row(p).copy();
        nc::NdArray<double> right = vertices.row(p).copy();
        left = left + 50.*(vertices.row(mn_idx) - vertices.row(p));
        right = right + 50.*(vertices.row(mx_idx) - vertices.row(p));
        return nc::vstack({left, right});
    }

    int complement(const vector<int>& vertex, int v1, int v2) {
        for (auto& e: vertex) {
            if (e != v1 and e != v2) return e;
        }
        return -1;
    }

    bool isIllegal(const pair<int, int>& vertexId) {
        auto it = edge_to_faces.find(edges[vertexId.first][vertexId.second]);
        if (it == edge_to_faces.end() or it->second.size() < 2) return false;
        auto faces = it->second;
        int f1 = *faces.begin();
        int f2 = *next(faces.begin());

        int r = complement(face_to_vertex[f1], vertexId.first, vertexId.second);
        int k = complement(face_to_vertex[f2], vertexId.first, vertexId.second);

        double a1, a2;
        {
            auto a = vertices.row(vertexId.first) - vertices.row(r);
            auto b = vertices.row(vertexId.second) - vertices.row(r);
            a1 = nc::dot(a, b)[0];
            a1 = a1 / (nc::norm(a)[0] * nc::norm(b)[0]);
        }
        {
            auto a = vertices.row(vertexId.first) - vertices.row(k);
            auto b = vertices.row(vertexId.second) - vertices.row(k);
            a2 = nc::dot(a, b)[0];
            a2 = a2 / (nc::norm(a)[0] * nc::norm(b)[0]);
        }

        auto ans = nc::arccos(a1) + nc::arccos(a2);
        return ans >= nc::constants::pi;
    }

    void legalizeEdge(int r, const pair<int, int>& edgesId, int face_id) {
        if (isIllegal(edgesId)) {
            auto faces = edge_to_faces[edges[edgesId.first][edgesId.second]];
            int f1 = *faces.begin();
            int f2 = *next(faces.begin());

            int r_ = complement(face_to_vertex[f1], edgesId.first, edgesId.second);
            int k = complement(face_to_vertex[f2], edgesId.first, edgesId.second);
            if (r == r_) r = r_;
            else if (r == k) k = r_;
            else {
                cerr << "error the vertex r is not in the cuadrilatero" << endl;
            }

            eraseTriangle(f1);
            eraseTriangle(f2);

            int face1 = addTriangle(r, k, edgesId.first);
            int face2 = addTriangle(r, k, edgesId.second);
            legalizeEdge(r, {k, edgesId.first}, face1);
            legalizeEdge(r, {k, edgesId.second}, face2);
        }
    }

    void addArise(int v1, int v2) {
        auto it = edges.find(v1);
        if (it == edges.end()) {
            edges[v1][v2] = edgeId;
            edges[v2][v1] = edgeId;
            edgeId++;
        }
        else {
            auto it2 = it->second.find(v2);
            if (it2 == it->second.end()) {
                edges[v1][v2] = edgeId;
                edges[v2][v1] = edgeId;
                edgeId++;
            }
        }
    }
    int addTriangle(int v1, int v2, int v3) {
        addArise(v1, v2);
        addArise(v2, v3);
        addArise(v3, v1);

        edge_to_faces[edges[v1][v2]].insert(faceId);
        edge_to_faces[edges[v2][v3]].insert(faceId);
        edge_to_faces[edges[v3][v1]].insert(faceId);

        face_to_vertex[faceId] = {v1, v2, v3};
        int tmp = faceId;
        faceId++;
        return tmp;
    }

    void eraseTriangle(int face) {
        auto vertexes = face_to_vertex[face];
        edge_to_faces[edges[vertexes[0]][vertexes[1]]].erase(face);
        edge_to_faces[edges[vertexes[1]][vertexes[2]]].erase(face);
        edge_to_faces[edges[vertexes[2]][vertexes[0]]].erase(face);

        face_to_vertex.erase(face);
    }

    pair<int, pair<int, int>> containInTriangle(int v){
        for (auto& [k, vertex]: face_to_vertex) {
            int l = vertex[0];
            int m = vertex[1];
            int n = vertex[2];
            auto res = this->contain(v, l, m, n);
            if (v == 11) cout << res.first << " " << res.second << endl;
            if (res.second != -1) return {k, res};
        }
        return {-1, {-1, -1}};
    }

    pair<int, int> contain(int v, int v1, int v2, int v3) {
        auto res1 = nc::cross(vertices.row(v2) - vertices.row(v1), vertices.row(v) - vertices.row(v1))[0];
        auto res2 = nc::cross(vertices.row(v3) - vertices.row(v2), vertices.row(v) - vertices.row(v2))[0];
        auto res3 = nc::cross(vertices.row(v1) - vertices.row(v3), vertices.row(v) - vertices.row(v3))[0];

        if ((res1 <= 0 && res2 <= 0 && res3 <= 0) or (res1 >= 0 && res2 >= 0 && res3 >= 0)){
            if (res1 == 0) return {v1, v2};
            if (res2 == 0) return {v2, v3};
            if (res3 == 0) return {v3, v1};
            return {-1, 0};
        }
        else return {-1, -1};

    }

    void print() {
        for (auto& [k, face]: face_to_vertex) {
            cout << to_string(k) << ":\t" << char('a' + face[0]) << " " << char('a' + face[1])  << " " << char('a' +
            face[2])
            << endl;
        }
    }
    void write() {
        fstream f("data.txt", ios::out);
        for (auto& [k, face]: face_to_vertex) {
            f << tt(vertices.row(face[0])) << " " << tt(vertices.row(face[1])) << " " << tt(vertices.row(face[2])) <<
            endl;
        }
    }
};



int main(int c, char *argv[]) {
    nc::NdArray<double> points;
    if (c > 1) {
        clog << "reading file " << argv[1] << endl;
        fstream f(argv[1], ios::in);
        vector<vector<double>> tmp;
        if (!f.is_open()) {
            cerr << "cant open" << endl;
        }
        while (true) {
            double a = 0, b = 0;
            f >> a >> b;
            clog << "readed " << a  <<  " " << b << endl;
            tmp.push_back({a, b});
            if (f.eof()) break;
        };
        points = nc::NdArray<double>(tmp);
    }
    points.shape().print();
    DelaunayTriangulation dt(points);
    dt.print();
    dt.write();
    return 0;
}
