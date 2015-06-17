#include <iostream>
#include <fstream>
#include <getopt.h>
#include "sdsl/bit_vectors.hpp"
#include "stream.hpp"
#include "cpp/vg.pb.h"

using namespace std;
using namespace sdsl;

class SuccinctGraph {
public:
    // build up interface here
    void add_node(Node& node);
    void add_edge(Edge& edge);
    Node node(int64_t rank); // id?
    Edge edge(int64_t rank);
    Path path(string& name);
    Graph neighborhood(int64_t rank, int32_t steps);
    Graph range(int64_t rank1, int64_t rank2);
    Graph region(string& path_name, int64_t start, int64_t stop);
private:
    // sequence/integer vector
    int_vector<> sequence;
    // node starts in sequence
    // rank_1(i) = node id
    // select_1(node_id) = i
    bit_vector node_in_seq;
    int_vector<> edges;
    bit_vector node_in_edges;
    map<string, bit_vector> paths;
};

void load_input_graph(string& input, list<& small_graph) {
    Person person1;
    person1.set_name("Ralph");
    person1.set_id(42);
    person1.set_email("post please");
    Person person2;
    person2.set_name("Sal");
    person2.set_id(99929999);
    person2.set_email("morse code");
    ofstream out;
    out.open("test.stream");
    vector<Person> people = { person1, person2 };
    stream::write_buffered(out, people, 0);
    out.close();
    ifstream in;
    in.open(input);
    function<void(Person&)> lambda = [](Person& p) {
        cout << p.name() << endl;
        cout << p.id() << endl;
    };
    stream::for_each(in, lambda);
    in.close();
}


int main(int argc, char** argv) {

    // what we need
    // process the graph to build up our entities
    // we need to get some statistics about the graph
    // how many nodes are there, how many edges, how much sequence
    // how big is the alphabet?
    // 
    // 
    // convert graph sequence into
    // sequence -> int_vector
    // edges into a list of targets and a list of sources
    // for each list we keep an identical-length bitvector which provides an id index over the edges
    
    bit_vector b(10000000, 0);
    b[8] = 1;
    rank_support_v<> rb(&b);

    cout<<rb(8)<<endl;
    cout<<rb(9)<<endl;

    cout<< "size of b in MB: " << size_in_mega_bytes(b)<< endl;
    cout<< "size of rb in MB: " << size_in_mega_bytes(rb)<< endl;

    rrr_vector<127> rrrb(b);
    rrr_vector<127>::rank_1_type rank_rrrb(&rrrb);
    cout<<rank_rrrb(8)<<endl;
    cout<<rank_rrrb(9)<<endl;

    cout<< "size of rrrb in MB: " << size_in_mega_bytes(rrrb)<< endl;
    cout<< "size of rank_rrrb in MB: " << size_in_mega_bytes(rank_rrrb)<< endl;


    rrr_vector<127>::select_1_type select_rrrb(&rrrb);
    cout<<"position of first one in b: "<<select_rrrb(1)<<endl;

    bit_vector x;
    util::assign(x, bit_vector(10000000,1));

    int_vector<> v(100, 5, 7);

    cout<<"v[5]="<<v[5]<<endl;
    v[5]=120;
    cout<<"v[5]="<<v[5]<<endl;


    int_vector<32> w(100, 4);

    write_structure<JSON_FORMAT>(rrrb, cout);
    cout<<endl;
}
