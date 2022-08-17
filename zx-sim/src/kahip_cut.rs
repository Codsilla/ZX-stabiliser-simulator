//Implementation of the extraction of a good vertex cut from a graph partion minimizing the number of edges

use std::collections::HashMap;


use quizx::{vec_graph::Graph, hash_graph::GraphLike};

pub fn kahip_cut_finder(g :&Graph,epsilon : f64) -> Vec<usize> {

    if g.tcount()==0 {
        return vec![];
    }

    let nb_spider = g.num_vertices();

    let vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((s, *i - 1)) }).collect::<HashMap<_,_>>();
    let reverse_vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((*i - 1, s)) }).collect::<HashMap<_,_>>();
    let edges = g.edges().map(|x| (vmapping[&x.0] ,vmapping[&x.1] )).collect::<Vec<_>>();

    let mut weights = vec![];

    for v in 0..nb_spider {
        if *g.phase(reverse_vmapping[&v]).denom() == 4 { weights.push(1); }
        else { weights.push(0);}
    }

    let graph = kahip::Graph::new(nb_spider, edges).with_vertex_weights(weights);



    let bad_index_separator = graph.node_separator(2, epsilon, kahip::Mode::Fast, 1234, false);
    let mut separator = vec![];

    for x in bad_index_separator{

        separator.push(reverse_vmapping[&x]);

    }


    separator
}
