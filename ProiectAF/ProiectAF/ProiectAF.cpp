#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <limits.h>

using namespace std;

// O structura ajutatoare pentru a retine un vector de muchii cu costuri
struct cost {
    int first;
    int second;
    int weight;
};

class Graf
{
    int tip; // Tipul grafului: 0 - Neorientat, 1 - Orientat
    int nr_noduri; // Numarul de noduri
    int nr_muchii; // Numarul de muchii
    vector <vector <int>> lista_vecini; // Graful reprezentat prin lista vecinilor fiecarui nod
    vector <cost> lista_muchii; // Graful reprezentat printr-un vector de muchii cu costuri
    vector <vector <pair <int, int>>> lista_costuri; // Graful reprezentat prin lista vecinilor cu costuri 

public:
    //  Definirea mai multor constructori in functie de felul in care
    // vom vrea sa fie reprezentat graful dat.
    Graf() {
        tip = 0;
        nr_noduri = 0;
        nr_muchii = 0;
    }
    Graf(int tip, int nr_noduri, int nr_muchii) {
        this->tip = tip;
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
    }
    Graf(int tip, int nr_noduri, int nr_muchii, vector <vector<int>> lista_vecini) {
        this->tip = tip;
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
        this->lista_vecini.resize(nr_noduri + 1);
        for (int i = 1; i <= nr_noduri; i++)
            for (int j = 0; j < lista_vecini[i].size(); j++)
                this->lista_vecini[i].push_back(lista_vecini[i][j]);
    }
    Graf(int tip, int nr_noduri, int nr_muchii, vector <cost> lista_muchii) {
        this->tip = tip;
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
        this->lista_muchii = lista_muchii;
    }
    Graf(int tip, int nr_noduri, int nr_muchii, vector <vector <pair <int, int>>> lista_costuri) {
        this->tip = tip;
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
        this->lista_costuri = lista_costuri;
    }
    ~Graf() { }

    //  Supraincarcarea operatorilor de citire si afisare
    // nefolositi in problemele rezolvare, dar utili
    // in cazul in care se vrea o citire/afisare mai explicita
    friend istream& operator>>(istream& in, Graf& g);
    friend ostream& operator<<(ostream& out, const Graf& g);

private: // O colectie de functii ajutatoare celor principale
    void dfs_muchiiIntoarcere(int i, vector <int>& viz, vector <int>& niv, vector<int>& niv_min, vector <vector<int>>& solutie);
    void dfs_tare_conex(int i, vector <int>& onStack, vector <int>& niv, vector <int>& niv_min, int& index, stack <int>& noduri, int& nr_comp, vector <vector<int>>& solutie);
    void dfs_sortare(int nod_curent, vector <int>& viz, vector <int>& solutie);
    void dfs_biconex(int i, int& index, vector <int>& niv, vector <int>& niv_min, vector <int>& tata, stack <pair<int, int>>& muchii, int& nr_comp, vector <vector <pair<int, int>>>& solutie);
    int reprez_kruskal(int nod, vector <int>& tata);
    void reuneste_kruskal(int ru, int rv, vector <int>& tata, vector <int>& h);
    int BFS_darb(int& start);
    bool bfs_flux(int s, int t, vector <int>& tata, vector <vector <int>>& flux);
    void euler(int nod, vector <int>& viz, vector <vector <pair <int, int>>>& lista_vec, vector <int>& solutie);
    void dfs_hamilton(int nod, vector <int>& viz, int nr_noduri_vizitate, int cost, int& cost_min);
    bool bfs_cuplaj(vector <int> &pereche_U, vector <int>& pereche_V, vector <int>& dist);
    bool dfs_cuplaj(int nod, vector <int>& pereche_U, vector <int>& pereche_V, vector <int>& dist);

public: // Functiile apelabile pentru un graf dat
    vector <int> BFS(int start);
    void DFS(int start, vector <int>& viz);
    int comp_conexe();
    bool havel_hakimi(int seq_len, vector <int> seq);
    vector <vector <int>> muchii_critice();
    vector <vector <int>> tare_conex();
    vector <int> sortare_top();
    vector <vector <pair <int, int>>> biconex();
    vector <pair <int, int>> kruskal();
    vector <bool> disjoint();
    vector <int> bellman_ford();
    vector <int> dijkstra();
    void floyd_warshall(vector <vector <int>>& matrice_costuri);
    int diametru();
    int edmond_karp(vector <vector <int>>& flux);
    vector <int> ciclueuler(vector <vector <pair <int, int>>>& lista_vec);
    int cicluhamilton();
    vector <pair <int, int>> hopcroft(int n, int m);
};

istream& operator>>(istream& in, Graf& g)
{
    cout << "\nIntroduceti tipul grafului (0-Neorientat 1-Orientat):";
    in >> g.tip;
    if (g.tip == 0) {
        cout << "\nIntroduceti numarul de noduri: ";
        cin >> g.nr_noduri;
        cout << "\nIntroduceti numarul de muchii: ";
        cin >> g.nr_muchii;
        cout << "\nIntroduceti lista muchiilor: ";
        g.lista_vecini.resize(g.nr_noduri + 1);
        for (int i = 0; i < g.nr_muchii; i++) {
            int x, y;
            cin >> x >> y;
            g.lista_vecini[x].push_back(y);
            g.lista_vecini[y].push_back(x);
        }
    }
    else if (g.tip == 1) {
        cout << "\nIntroduceti numarul de noduri: ";
        cin >> g.nr_noduri;
        cout << "\nIntroduceti numarul de arce: ";
        cin >> g.nr_muchii;
        cout << "\nIntroduceti lista arcelor: ";
        g.lista_vecini.resize(g.nr_noduri + 1);
        for (int i = 0; i < g.nr_muchii; i++) {
            int x, y;
            cin >> x >> y;
            g.lista_vecini[x].push_back(y);
        }
    }
    else {
        cout << "\nTip gresit";
    }

    return in;
}
ostream& operator<<(ostream& out, const Graf& g)
{
    cout << "\nTip graf: ";
    if (g.tip == 0) {
        cout << "neorientat";
        cout << "\nNumar de noduri: " << g.nr_noduri;
        cout << "\nNumar de muchii: " << g.nr_muchii;
        cout << "\nLista de vecini: ";
        for (int i = 1; i <= g.nr_noduri; i++)
        {
            for (int j = 0; j < g.lista_vecini[i].size(); j++)
                cout << g.lista_vecini[i][j] << " ";
            cout << "\n";
        }
    }
    else if (g.tip == 1) {
        cout << "orientat";
        cout << "\nNumar de noduri: " << g.nr_noduri;
        cout << "\nNumar de arce: " << g.nr_muchii;
        cout << "\nLista de vecini: ";
        for (int i = 1; i <= g.nr_noduri; i++)
        {
            for (int j = 0; j < g.lista_vecini[i].size(); j++)
                cout << g.lista_vecini[i][j] << " ";
            cout << "\n";
        }
    }

    return out;
}

// O colectie de metode care rezolva problemele date
void infoarena_bfs() {
    //  Se da un graf orientat si un nod de start.
    //  Se cere sa se afiseze distanta de la nodul de start pana la
    // celelate noduri sau -1 in cazul in care nu se poate ajunge la el.

    ifstream in("bfs.in");
    ofstream out("bfs.out");

    int n, m, s;
    vector <vector <int>> v;

    in >> n >> m >> s;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
    }
    Graf g(1, n, m, v);

    vector <int> solutie;
    solutie = g.BFS(s);
    for (int i = 1; i <= n; i++)
        out << solutie[i] << " ";

    in.close();
    out.close();
}
void infoarena_dfs() {
    //  Se da un graf neorientat.
    //  Se cere sa se afiseze numarul de componente conexe ale grafului
    // folosindu-se un DFS.

    ifstream in("dfs.in");
    ofstream out("dfs.out");

    int n, m;
    vector <vector <int>> v;

    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
        v[y].push_back(x);
    }
    Graf g(0, n, m, v);

    out << g.comp_conexe();

    in.close();
    out.close();
}
void rezolva_havel() {
    //  Se da o secventa de numere.
    //  Se cere sa se afiseze daca se poate sau nu sa se formeze
    // un graf folosind secventa de numere data.

    ifstream in("havel.in");
    ofstream out("havel.out");

    int n;
    vector <int> v;

    in >> n;
    for (int i = 0; i < n; i++) {
        int x;
        in >> x;
        v.push_back(x);
    }
    Graf g;

    if (g.havel_hakimi(n, v) == true)
        out << "Da";
    else
        out << "Nu";

    in.close();
    out.close();
}
void leetcode_criticalCon() {
    //  Se da un graf neorientat.
    //  Se cere sa se afiseze numarul de muchii critice alea grafului
    // si lista acestora.

    ifstream in("critcon.in");
    ofstream out("critcon.out");

    int n, m;
    vector <vector <int>> v;

    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
        v[y].push_back(x);
    }
    Graf g(0, n, m, v);

    vector <vector <int>> solutie;
    solutie = g.muchii_critice();
    for (int i = 0; i < solutie.size(); i++)
        out << solutie[i][0] << " " << solutie[i][1] << "\n";

    in.close();
    out.close();
}
void infoarena_ctc() {
    //  Se da un graf orientat.
    //  Se cere sa se afiseze numarul de componente tare conexe ale grafului,
    // precum si care sunt nodurile care le formeaza.

    ifstream in("ctc.in");
    ofstream out("ctc.out");

    int n, m;
    vector <vector <int>> v;

    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
    }
    Graf g(1, n, m, v);

    vector <vector <int>> solutie;
    solutie = g.tare_conex();
    int nr_comp = solutie[0][0];

    out << nr_comp << "\n";
    for (int i = 1; i <= nr_comp; i++) {
        for (int j = 0; j < solutie[i].size(); j++)
            out << solutie[i][j] << " ";
        out << "\n";
    }

    in.close();
    out.close();
}
void infoarena_sortaret() {
    //  Se da un graf orientat.
    //  Se cere sa se afiseze sortarea topologica a nodurilor grafului dat.

    ifstream in("sortaret.in");
    ofstream out("sortaret.out");

    int n, m;
    vector <vector <int>> v;

    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
    }
    Graf g(1, n, m, v);

    vector <int> solutie;
    solutie = g.sortare_top();
    for (int i = solutie.size() - 1; i >= 0; i--)
        out << solutie[i] << " ";

    in.close();
    out.close();
}
void infoarena_biconex() {
    //  Se da un graf neorientat.
    //  Se cere sa se afiseze numarul de componente biconexe ale grafului,
    // precum si nodurile care compun fiecare componenta.

    ifstream in("biconex.in");
    ofstream out("biconex.out");

    int n, m;
    vector <vector <int>> v;

    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
        v[y].push_back(x);
    }
    Graf g(0, n, m, v);

    vector <vector <pair <int, int>>> solutie;
    solutie = g.biconex();
    int nr_comp = solutie[0][0].first;
    out << nr_comp << "\n";
    for (int i = 1; i <= nr_comp; i++) {
        vector <int> noduri(n + 1, 0);

        for (int j = 0; j < solutie[i].size(); j++) {
            if (noduri[solutie[i][j].first] == 0) {
                out << solutie[i][j].first << " ";
                noduri[solutie[i][j].first] = 1;
            }
            if (noduri[solutie[i][j].second] == 0) {
                out << solutie[i][j].second << " ";
                noduri[solutie[i][j].second] = 1;
            }
        }
        out << "\n";
    }

    in.close();
    out.close();
}
void infoarena_apm() {
    //  Se da un graf neorientat cu costuri.
    //  Se cere sa se afiseze costul total al arborelui partial de cost minim,
    // numarul de muchii ale acestuia, precum si vectorul muchiilor ce fac parte
    // din el.

    ifstream in("apm.in");
    ofstream out("apm.out");

    int n, m;
    vector <cost> v;

    in >> n >> m;
    for (int i = 0; i < m; i++) {
        int x, y, c;
        in >> x >> y >> c;
        v.push_back({ x,y,c });
    }
    Graf g(0, n, m, v);

    vector <pair <int, int>> solutie;
    solutie = g.kruskal();
    int cost_total = solutie[0].first;
    int nr_muchiiSol = solutie[0].second;
    out << cost_total << "\n" << nr_muchiiSol << "\n";
    for (int i = 1; i <= nr_muchiiSol; i++)
        out << solutie[i].first << " " << solutie[i].second << "\n";

    in.close();
    out.close();
}
void infoarena_disjoint() {
    //  Se dau mai multe triplete de numere de format (cod, x, y) cu semnificatia urmatoare:
    // cod = 1  =>  Se cere sa se reuneasca multimile in care se afla x si y.
    // cod = 2  =>  Se cere sa se afiseze daca x si y fac parte din aceiasi multime sau nu.

    ifstream in("disjoint.in");
    ofstream out("disjoint.out");

    int n, m;
    vector <cost> v;

    in >> n >> m;
    for (int i = 0; i < m; i++) {
        int x, y, cod;
        in >> cod >> x >> y;
        v.push_back({ x,y,cod });
    }

    //  Pentru aceasta problema am folosit atributele clasei cu urmatoarea semnificatie:
    // nr_noduri = numarul de multimi
    // nr_muchii = numarul de operatii
    // lista_muchii = lista tripletelor de forma (cod, x, y)
    Graf g(0, n, m, v);

    vector <bool> solutie;
    solutie = g.disjoint();
    for (int i = 0; i < solutie.size(); i++)
        if (solutie[i] == true)
            out << "DA\n";
        else
            out << "NU\n";

    in.close();
    out.close();
}
void infoarena_bellmanford() {
    //  Se da un graf orientat cu costuri, care pot fi si negative.
    //  Se cere sa se afiseze costul minim al unui lant de la nodul 1 la fiecare nod,
    // sau un mesaj in cazul in care graful are cicluri negative.

    ifstream in("bellmanford.in");
    ofstream out("bellmanford.out");

    int n, m;
    vector <vector <pair<int, int>>> v;

    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y, c;
        in >> x >> y >> c;
        v[x].push_back(make_pair(y, c));
    }
    Graf g(1, n, m, v);

    vector <int> solutie;
    solutie = g.bellman_ford();
    if (!solutie.empty())
        for (int i = 2; i <= n; i++)
            out << solutie[i] << " ";
    else
        out << "Ciclu negativ!";

    in.close();
    out.close();
}
void infoarena_dijkstra() {
    //  Se da un graf orientat cu costuri pozitive.
    //  Se cere sa se afiseze costul minim al unui lant de la nodul 1 la fiecare nod.

    ifstream in("dijkstra.in");
    ofstream out("dijkstra.out");

    int n, m;
    vector <vector <pair <int, int>>> v;

    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y, c;
        in >> x >> y >> c;
        v[x].push_back(make_pair(y, c));
    }
    Graf g(1, n, m, v);

    vector <int> solutie;
    solutie = g.dijkstra();
    int inf = INT_MAX;
    for (int i = 2; i <= n; i++)
        if (solutie[i] == inf)
            out << 0 << " ";
        else
            out << solutie[i] << " ";

    in.close();
    out.close();
}
void infoarena_floydwarshall() {
    //  Se da un graf orientat cu costuri.
    //  Se cere sa se afiseze matricea drumurilor de cost minim.

    ifstream in("royfloyd.in");
    ofstream out("royfloyd.out");

    int n;
    in >> n;
    vector <vector <int>> d(n + 1, vector <int>(n + 1, 0));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            in >> d[i][j];
    Graf g(1, n, -1);

    g.floyd_warshall(d);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++)
            out << d[i][j] << " ";
        out << "\n";
    }

    in.close();
    out.close();
}
void infoarena_diametru() {
    //  Se da un arbore si se cere afisarea diametrului sau.

    ifstream in("darb.in");
    ofstream out("darb.out");

    int n, m;
    in >> n;
    m = n - 1;
    vector <vector <int>> v;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
        v[y].push_back(x);
    }
    Graf g(0, n, m, v);

    int d = g.diametru();
    out << d;

    in.close();
    out.close();
}
void infoarena_maxflow() {
    //  Se da un graf orientat cu fluxuri.
    //  Se cere afisarea fluxului maxim ce poate sa ajunga din nodul 1 la nodul N.

    ifstream in("maxflow.in");
    ofstream out("maxflow.out");

    int n, m;
    vector <vector <int>> v;
    in >> n >> m;
    v.resize(n + 1);
    vector <vector <int>> c(n + 1, vector <int>(n + 1, 0));
    for (int i = 0; i < m; i++) {
        int x, y, z;
        in >> x >> y >> z;
        v[x].push_back(y);
        v[y].push_back(x);
        c[x][y] = z;
    }
    Graf g(1, n, m, v);

    out << g.edmond_karp(c);

    in.close();
    out.close();
}
void infoarena_euler() {
    //  Se da un multigraf.
    //  Se cere sa se afiseze nodurile care alcatuiesc un ciclu eulerian.

    ifstream in("ciclueuler.in");
    ofstream out("ciclueuler.out");

    int n, m;
    vector <vector <pair <int, int>>> v;
    in >> n >> m;
    v.resize(n + 1);
    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(make_pair(y, i + 1));
        v[y].push_back(make_pair(x, i + 1));
    }
    Graf g(0, n, m);

    vector <int> sol = g.ciclueuler(v);
    for (int i = 0; i < sol.size(); i++)
        out << sol[i] << " ";

    in.close();
    out.close();
}
void infoarena_hamilton() {
    //  Se da un graf orientat cu costuri.
    //  Se cere sa se afiseze costul minim al unui ciclu hamiltonian.

    ifstream in("hamilton.in");
    ofstream out("hamilton.out");

    int n, m;
    vector <vector <pair <int, int>>> v;
    in >> n >> m;
    v.resize(n);
    for (int i = 0; i < m; i++) {
        int x, y, c;
        in >> x >> y >> c;
        v[x].push_back(make_pair(y, c));
    }
    Graf g(1, n, m, v);

    int sol = g.cicluhamilton();
    if (sol == 0)
        out << "Nu exista solutie";
    else
        out << sol;

    in.close();
    out.close();
}
void infoarena_bipartit() {
    //  Se da un graf neorientat bipartit.
    //  Se cere sa se afiseze un cuplaj maxim din graf.

    ifstream in("cuplaj.in");
    ofstream out("cuplaj.out");

    int n, m, e;
    in >> n >> m >> e;
    vector <vector <int>> v;
    v.resize(n + 1);
    for (int i = 0; i < e; i++) {
        int x, y;
        in >> x >> y;
        v[x].push_back(y);
    }
    Graf g(0, n, e, v);

    vector <pair <int, int>> sol;
    sol = g.hopcroft(n, m);

    out << sol.back().first << "\n";
    for (int i = 0; i < sol.size() - 1; i++)
        out << sol[i].first << " " << sol[i].second << "\n";

    in.close();
    out.close();
}

int main()
{
    //  Pentru rezolvarea unei anumite probleme date pe infoarena/leetcode
    // doar apelati metoda respectiva.
    //  Datele de intrare si cele de iesire se vor scrie in fisierele 
    // corespunzatoare problemei.
    infoarena_hamilton();
    return 0;
}



vector <int> Graf::BFS(int start) {
    queue <int> bfs_tree;
    vector <int> viz(nr_noduri + 1, -1);

    bfs_tree.push(start);
    viz[start] = 0;

    while (!bfs_tree.empty()) {
        int nod_curent = bfs_tree.front();
        bfs_tree.pop();

        for (int i = 0; i < lista_vecini[nod_curent].size(); i++)
            if (viz[lista_vecini[nod_curent][i]] == -1) {
                viz[lista_vecini[nod_curent][i]] = viz[nod_curent] + 1;
                bfs_tree.push(lista_vecini[nod_curent][i]);
            }
    }

    return viz;
}
void Graf::DFS(int start, vector <int>& viz) {
    stack <int> dfs_tree;

    dfs_tree.push(start);
    viz[start] = 1;

    while (!dfs_tree.empty()) {
        int nod_curent = dfs_tree.top();
        int vecin = 0;
        for (int i = 0; i < lista_vecini[nod_curent].size(); i++)
            if (viz[lista_vecini[nod_curent][i]] == 0) {
                vecin = lista_vecini[nod_curent][i];
                break;
            }
        if (vecin == 0)
            dfs_tree.pop();
        else {
            dfs_tree.push(vecin);
            viz[vecin] = 1;
        }
    }

}
int Graf::comp_conexe() {
    vector <int> viz(nr_noduri + 1, 0);
    int nr_comp = 0;

    for (int i = 1; i <= this->nr_noduri; i++)
        if (viz[i] == 0) {
            DFS(i, viz);
            nr_comp++;
        }

    return nr_comp;
}
bool Graf::havel_hakimi(int seq_len, vector <int> seq) {
    int sum = 0;
    for (int i = 0; i < seq_len; i++) {
        sum = sum + seq[i];
        if (seq[i] > seq_len - 1)
            return false;
    }
    if (sum % 2 != 0)
        return false;
    else {
        sort(seq.begin(), seq.end(), greater<int>());
        while (seq[0] != 0) {
            int nod_curent = seq[0];
            seq.erase(seq.begin());
            for (int i = 0; i < nod_curent; i++) {
                seq[i]--;
                if (seq[i] < 0)
                    return false;
            }
            sort(seq.begin(), seq.end(), greater<int>());
        }
        return true;
    }
}
void Graf::dfs_muchiiIntoarcere(int i, vector <int>& viz, vector <int>& niv, vector<int>& niv_min, vector <vector<int>>& muchii_crit) {
    viz[i] = 1;
    niv_min[i] = niv[i];
    for (int j = 0; j < lista_vecini[i].size(); j++)
        if (viz[lista_vecini[i][j]] == 0) {
            niv[lista_vecini[i][j]] = niv[i] + 1;
            dfs_muchiiIntoarcere(lista_vecini[i][j], viz, niv, niv_min, muchii_crit);
            niv_min[i] = min(niv_min[i], niv_min[lista_vecini[i][j]]);
            if (niv_min[lista_vecini[i][j]] > niv[i]) {
                muchii_crit.push_back({ i,lista_vecini[i][j] });
            }
        }
        else {
            if (niv[lista_vecini[i][j]] < niv[i] - 1)
                niv_min[i] = min(niv_min[i], niv[lista_vecini[i][j]]);
        }
}
vector <vector <int>> Graf::muchii_critice() {
    vector <vector<int>> muchii_crit;
    vector <int> viz(nr_noduri + 1, 0);
    vector <int> niv;
    vector <int> niv_min;

    niv.resize(nr_noduri + 1);
    niv_min.resize(nr_noduri + 1);
    niv[1] = 1;
    dfs_muchiiIntoarcere(1, viz, niv, niv_min, muchii_crit);

    return muchii_crit;
}
void Graf::dfs_tare_conex(int i, vector <int>& onStack, vector <int>& niv, vector <int>& niv_min, int& index, stack <int>& noduri, int& nr_comp, vector <vector<int>>& solutie) {
    niv[i] = index;
    niv_min[i] = index;
    index++;
    noduri.push(i);
    onStack[i] = 1;

    for (int j = 0; j < lista_vecini[i].size(); j++)
        if (niv[lista_vecini[i][j]] == 0) {
            dfs_tare_conex(lista_vecini[i][j], onStack, niv, niv_min, index, noduri, nr_comp, solutie);
            niv_min[i] = min(niv_min[i], niv_min[lista_vecini[i][j]]);
        }
        else if (onStack[lista_vecini[i][j]] == 1)
            niv_min[i] = min(niv_min[i], niv[lista_vecini[i][j]]);

    if (niv_min[i] == niv[i]) {
        nr_comp++;
        while (noduri.top() != i) {
            solutie[nr_comp].push_back(noduri.top());
            onStack[noduri.top()] = 0;
            noduri.pop();
        }
        solutie[nr_comp].push_back(noduri.top());
        onStack[noduri.top()] = 0;
        noduri.pop();
    }
}
vector <vector <int>> Graf::tare_conex() {
    vector <int> niv(nr_noduri + 1, 0);
    vector <int> niv_min;
    stack <int> noduri;
    vector <int> onStack(nr_noduri + 1, 0);
    vector <vector<int>> solutie;

    int nr_comp = 0;
    niv_min.resize(nr_noduri + 1);
    solutie.resize(nr_noduri + 1);
    int index = 1;
    for (int i = 1; i <= nr_noduri; i++)
        if (niv[i] == 0)
            dfs_tare_conex(i, onStack, niv, niv_min, index, noduri, nr_comp, solutie);

    solutie[0].push_back(nr_comp);
    return solutie;
}
void Graf::dfs_sortare(int nod_curent, vector <int>& viz, vector <int>& solutie) {
    viz[nod_curent] = 1;
    for (int i = 0; i < lista_vecini[nod_curent].size(); i++)
        if (viz[lista_vecini[nod_curent][i]] == 0)
            dfs_sortare(lista_vecini[nod_curent][i], viz, solutie);
    solutie.push_back(nod_curent);
}
vector <int> Graf::sortare_top() {
    vector <int> solutie;
    vector <int> viz(nr_noduri + 1, 0);

    for (int i = 1; i <= nr_noduri; i++)
        if (viz[i] == 0)
            dfs_sortare(i, viz, solutie);

    return solutie;
}
void Graf::dfs_biconex(int i, int& index, vector <int>& niv, vector <int>& niv_min, vector <int>& tata, stack <pair<int, int>>& muchii, int& nr_comp, vector <vector <pair<int, int>>>& solutie) {
    niv[i] = index;
    niv_min[i] = index;
    index++;
    int nr_copii = 0;

    for (int j = 0; j < lista_vecini[i].size(); j++) {
        if (niv[lista_vecini[i][j]] == 0) {
            nr_copii++;
            tata[lista_vecini[i][j]] = i;
            muchii.push(make_pair(i, lista_vecini[i][j]));

            dfs_biconex(lista_vecini[i][j], index, niv, niv_min, tata, muchii, nr_comp, solutie);

            niv_min[i] = min(niv_min[i], niv_min[lista_vecini[i][j]]);
            if ((niv[i] == 1 && nr_copii > 1) || (niv[i] > 1 && niv_min[lista_vecini[i][j]] >= niv[i])) {
                nr_comp++;
                while (muchii.top().first != i || muchii.top().second != lista_vecini[i][j]) {
                    solutie[nr_comp].push_back(muchii.top());
                    muchii.pop();
                }
                solutie[nr_comp].push_back(muchii.top());
                muchii.pop();
            }
        }
        else if (lista_vecini[i][j] != tata[i]) {
            niv_min[i] = min(niv_min[i], niv[lista_vecini[i][j]]);
            if (niv[lista_vecini[i][j]] < niv[i])
                muchii.push(make_pair(i, lista_vecini[i][j]));
        }
    }
}
vector <vector <pair <int, int>>> Graf::biconex() {
    vector <int> niv(nr_noduri + 1, 0);
    vector <int> niv_min;
    vector <int> tata(nr_noduri + 1, 0);
    stack <pair <int, int>> muchii;
    vector <vector <pair <int, int>>> solutie;
    niv_min.resize(nr_noduri + 1);
    solutie.resize(nr_noduri);

    int index = 1;
    int nr_comp = 0;

    for (int i = 1; i <= nr_noduri; i++)
        if (niv[i] == 0)
            dfs_biconex(i, index, niv, niv_min, tata, muchii, nr_comp, solutie);

    if (!muchii.empty()) {
        nr_comp++;
        while (!muchii.empty()) {
            solutie[nr_comp].push_back(muchii.top());
            muchii.pop();
        }
    }
    solutie[0].push_back(make_pair(nr_comp, 0));

    return solutie;
}
int Graf::reprez_kruskal(int nod, vector <int>& tata) {
    while (tata[nod] != 0)
        nod = tata[nod];
    return nod;
}
void Graf::reuneste_kruskal(int ru, int rv, vector <int>& tata, vector <int>& h) {
    if (h[ru] > h[rv])
        tata[rv] = ru;
    else {
        tata[ru] = rv;
        if (h[ru] == h[rv])
            h[rv]++;
    }
}
bool criteriuSort(cost c1, cost c2) { return (c1.weight < c2.weight); }
vector <pair <int, int>> Graf::kruskal() {
    vector <pair <int, int>> solutie;
    solutie.resize(nr_noduri);
    vector <int> tata(nr_noduri + 1, 0);
    vector <int> h(nr_noduri + 1, 0);

    sort(lista_muchii.begin(), lista_muchii.end(), criteriuSort);

    int cost_total = 0;
    int nr_muchiiSol = 0;
    for (int i = 0; (i < lista_muchii.size()) && (nr_muchiiSol < nr_noduri - 1); i++) {
        int ru = reprez_kruskal(lista_muchii[i].first, tata);
        int rv = reprez_kruskal(lista_muchii[i].second, tata);
        if (ru != rv) {
            nr_muchiiSol++;
            solutie[nr_muchiiSol].first = lista_muchii[i].first;
            solutie[nr_muchiiSol].second = lista_muchii[i].second;
            cost_total = cost_total + lista_muchii[i].weight;
            reuneste_kruskal(ru, rv, tata, h);
        }
    }

    solutie[0].first = cost_total;
    solutie[0].second = nr_muchiiSol;
    return solutie;
}
vector <bool> Graf::disjoint() {
    vector <bool> solutie;
    vector <int> tata(nr_noduri + 1, 0);
    vector <int> h(nr_noduri + 1, 0);

    for (int i = 0; i < nr_muchii; i++)
        if (lista_muchii[i].weight == 1) {
            int ru = reprez_kruskal(lista_muchii[i].first, tata);
            int rv = reprez_kruskal(lista_muchii[i].second, tata);
            reuneste_kruskal(ru, rv, tata, h);
        }
        else if (lista_muchii[i].weight == 2) {
            int ru = reprez_kruskal(lista_muchii[i].first, tata);
            int rv = reprez_kruskal(lista_muchii[i].second, tata);
            if (ru == rv)
                solutie.push_back(true);
            else
                solutie.push_back(false);
        }

    return solutie;
}
vector <int> Graf::bellman_ford() {
    int inf = INT_MAX;
    vector <int> distanta(nr_noduri + 1, inf);
    queue <int> coada;
    vector <bool> inCoada(nr_noduri + 1, false);
    vector <int> pasi(nr_noduri + 1, 0);

    distanta[1] = 0;
    coada.push(1);
    inCoada[1] = true;
    while (!coada.empty()) {
        int u = coada.front();
        coada.pop();
        inCoada[u] = false;
        for (int i = 0; i < lista_costuri[u].size(); i++) {
            int v = lista_costuri[u][i].first;
            int c = lista_costuri[u][i].second;
            if (distanta[u] + c < distanta[v]) {
                distanta[v] = distanta[u] + c;
                pasi[v]++;
                if (pasi[v] == nr_noduri)
                    return {};
                if (inCoada[v] == false) {
                    coada.push(v);
                    inCoada[v] = true;
                }
            }
        }
    }

    return distanta;
}
vector <int> Graf::dijkstra() {
    int inf = INT_MAX;
    vector <int> distanta(nr_noduri + 1, inf);
    priority_queue <pair <int, int>> min_heap;

    distanta[1] = 0;
    min_heap.push(make_pair(0, 1));
    while (!min_heap.empty()) {
        int u = min_heap.top().second;
        min_heap.pop();
        for (int i = 0; i < lista_costuri[u].size(); i++) {
            int v = lista_costuri[u][i].first;
            int c = lista_costuri[u][i].second;
            if (distanta[u] + c < distanta[v]) {
                distanta[v] = distanta[u] + c;
                min_heap.push(make_pair(-distanta[v], v));
            }
        }
    }

    return distanta;
}
void Graf::floyd_warshall(vector <vector <int>>& matrice_costuri) {
    for (int k = 1; k <= nr_noduri; k++)
        for (int i = 1; i <= nr_noduri; i++)
            for (int j = 1; j <= nr_noduri; j++)
                if (matrice_costuri[i][k] != 0 && matrice_costuri[k][j] != 0 && i != j)
                    if (matrice_costuri[i][j] > matrice_costuri[i][k] + matrice_costuri[k][j] || matrice_costuri[i][j] == 0)
                        matrice_costuri[i][j] = matrice_costuri[i][k] + matrice_costuri[k][j];
}
int Graf::BFS_darb(int& start) {
    vector <int> dist(nr_noduri + 1, 0);
    queue <int> coada;
    int diam;

    coada.push(start);
    dist[start] = 1;
    diam = dist[start];

    while (!coada.empty()) {
        int nod_curent = coada.front();
        start = nod_curent;
        coada.pop();
        for (int i = 0; i < lista_vecini[nod_curent].size(); i++)
            if (dist[lista_vecini[nod_curent][i]] == 0) {
                coada.push(lista_vecini[nod_curent][i]);
                dist[lista_vecini[nod_curent][i]] = dist[nod_curent] + 1;
                diam = dist[lista_vecini[nod_curent][i]];
            }
    }

    return diam;
}
int Graf::diametru() {
    int start = 1;
    int diam = BFS_darb(start);
    diam = BFS_darb(start);
    return diam;
}
bool Graf::bfs_flux(int s, int t, vector <int>& tata, vector <vector <int>>& flux) {
    vector <int> viz(nr_noduri + 1, 0);
    queue <int> coada;

    coada.push(s);
    viz[s] = 1;
    while (!coada.empty()) {
        int nod_curent = coada.front();
        coada.pop();
        for (int i = 0; i < lista_vecini[nod_curent].size(); i++) {
            if (viz[lista_vecini[nod_curent][i]] == 0 && flux[nod_curent][lista_vecini[nod_curent][i]] > 0) {
                if (lista_vecini[nod_curent][i] == t) {
                    tata[t] = nod_curent;
                    return true;
                }
                coada.push(lista_vecini[nod_curent][i]);
                viz[lista_vecini[nod_curent][i]] = 1;
                tata[lista_vecini[nod_curent][i]] = nod_curent;
            }
        }
    }
    return false;
}
int Graf::edmond_karp(vector <vector <int>>& flux) {
    int flux_max = 0;
    vector <int> tata(nr_noduri + 1, 0);
    while (bfs_flux(1, nr_noduri, tata, flux) == true) {
        int nod = nr_noduri;
        int flux_min = 0;
        while (nod != 1) {
            if (flux_min == 0 || flux[tata[nod]][nod] < flux_min)
                flux_min = flux[tata[nod]][nod];
            nod = tata[nod];
        }
        nod = nr_noduri;
        while (nod != 1) {
            flux[tata[nod]][nod] = flux[tata[nod]][nod] - flux_min;
            flux[nod][tata[nod]] = flux_min;
            nod = tata[nod];
        }
        flux_max = flux_max + flux_min;
    }
    return flux_max;
}
void Graf::euler(int nod, vector <int>& viz, vector <vector <pair <int, int>>>& lista_vec, vector <int>& solutie) {
    while (!lista_vec[nod].empty()) {
        int vecin = lista_vec[nod][lista_vec[nod].size() - 1].first;
        int nr_muchie = lista_vec[nod][lista_vec[nod].size() - 1].second;
        lista_vec[nod].pop_back();
        if (viz[nr_muchie] == 0) {
            viz[nr_muchie] = 1;
            euler(vecin, viz, lista_vec, solutie);
        }
    }
    solutie.push_back(nod);
}
vector <int> Graf::ciclueuler(vector <vector <pair <int, int>>>& lista_vec) {
    vector <int> solutie;

    for (int i = 1; i <= nr_noduri; i++)
        if (lista_vec[i].size() % 2 != 0)
            return { -1 };

    vector <int> viz(nr_muchii + 1, 0);
    euler(1, viz, lista_vec, solutie);

    return solutie;
}
void Graf::dfs_hamilton(int nod, vector <int>& viz, int nr_noduri_vizitate, int cost, int& cost_min) {
    viz[nod] = 1;
    for (int i = 0; i < lista_costuri[nod].size(); i++) {
        if (viz[lista_costuri[nod][i].first] == 0)
            dfs_hamilton(lista_costuri[nod][i].first, viz, nr_noduri_vizitate + 1, cost + lista_costuri[nod][i].second, cost_min);
        if (nr_noduri_vizitate == nr_noduri - 1 && lista_costuri[nod][i].first == 0)
            if (cost + lista_costuri[nod][i].second < cost_min || cost_min == 0)
                cost_min = cost + lista_costuri[nod][i].second;
    }
    viz[nod] = 0;
}
int Graf::cicluhamilton() {
    vector <int> viz(nr_noduri, 0);
    int cost_min = 0;
    dfs_hamilton(0, viz, 0, 0, cost_min);
    return cost_min;
}
bool Graf::bfs_cuplaj(vector <int>& pereche_U, vector <int>& pereche_V, vector <int>& dist) {
    queue <int> Q;
    int INF = INT_MAX;
    for (int i = 1; i < pereche_U.size(); i++)
        if (pereche_U[i] == 0) {
            dist[i] = 0;
            Q.push(i);
        }
        else
            dist[i] = INF;
    dist[0] = INF;
    while (!Q.empty()) {
        int nod = Q.front();
        Q.pop();
        if (dist[nod] < dist[0])
            for (int i = 0; i < lista_vecini[nod].size(); i++)
                if (dist[pereche_V[lista_vecini[nod][i]]] == INF) {
                    dist[pereche_V[lista_vecini[nod][i]]] = dist[nod] + 1;
                    Q.push(pereche_V[lista_vecini[nod][i]]);
                }
    }
    if (dist[0] != INF)
        return true;
    else
        return false;
}
bool Graf::dfs_cuplaj(int nod, vector <int>& pereche_U, vector <int>& pereche_V, vector <int>& dist) {
    int INF = INT_MAX;
    if (nod != 0) {
        for (int i = 0; i < lista_vecini[nod].size(); i++)
            if (dist[pereche_V[lista_vecini[nod][i]]] == dist[nod] + 1)
                if (dfs_cuplaj(pereche_V[lista_vecini[nod][i]], pereche_U, pereche_V, dist) == true) {
                    pereche_V[lista_vecini[nod][i]] = nod;
                    pereche_U[nod] = lista_vecini[nod][i];
                    return true;
                }
        dist[nod] = INF;
        return false;
    }
    return true;
}
vector <pair <int, int>> Graf::hopcroft(int nr_noduri_U, int nr_noduri_V) {
    vector <int> pereche_U(nr_noduri_U + 1, 0);
    vector <int> pereche_V(nr_noduri_V + 1, 0);
    vector <int> dist(nr_noduri_U + 1, 0);
    vector <pair <int, int>> sol;
    int nr_cmax = 0;
    while (bfs_cuplaj(pereche_U, pereche_V, dist) == true)
        for (int i = 1; i <= nr_noduri_U; i++)
            if (pereche_U[i] == 0 && dfs_cuplaj(i, pereche_U, pereche_V, dist) == true)
                nr_cmax++;
    for (int i = 1; i <= nr_noduri_U; i++)
        if (pereche_U[i] != 0)
            sol.push_back(make_pair(i, pereche_U[i]));
    sol.push_back(make_pair(nr_cmax, 0));
    return sol;
}