



#![allow(non_snake_case)]
use quizx::{vec_graph::Graph, hash_graph::GraphLike};
use kahypar::{self};
use std::collections::HashMap;

pub fn kaHyPar_cut_01(g :&Graph)-> Vec<usize> {
    kaHyPar_cut_finder(&g,0.1)
} 

pub fn kaHyPar_cut_finder(g :&Graph,epsilon : f64) -> Vec<usize> {

    if g.tcount()==0 {
        return vec![];
    }

    let nb_spider = g.num_vertices();
    let mut hyperedges = vec![Vec::<u32>::new();nb_spider];
    
    let vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((s, *i - 1)) }).collect::<HashMap<_,_>>();
    let reverse_vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((*i - 1, s)) }).collect::<HashMap<_,_>>();
    //let emapping = g.edges().scan(0, |i, (s, e, _)| { *i += 1; Some(((s, e), *i - 1)) }).collect::<HashMap<_,_>>();

    let edge_list = g.edges().map(|x| (x.0 ,x.1 )).collect::<Vec<_>>();

    for (i,(s,d)) in edge_list.iter().enumerate() {
        
        hyperedges[vmapping[s]].push(i as u32);
        hyperedges[vmapping[d]].push(i as u32);
    }


    let mut gh = kahypar::Hypergraph::from_edges(2,edge_list.len(),&hyperedges).build();

    let mut context = kahypar::Context::from_config("src/config.ini");
    let (_,part) = gh.partition(&mut context, epsilon);


    //build a list of all index reduced adjacent part for each vertex (with repetition)
    let mut adjacent_part = vec![Vec::<u32>::new();nb_spider];

    //TODO: pairs of bool instead
    for (i,(s,d)) in edge_list.iter().enumerate() {        
        adjacent_part[vmapping[s]].push(part[i] as u32);
        adjacent_part[vmapping[d]].push(part[i] as u32);
    }

    let reduced_index_cut = adjacent_part.iter().map(|x| is_all_same(x));


    //extract the cut with the right index
    let mut cut = Vec::new();
    for (ri,b) in reduced_index_cut.enumerate(){
        if !b {
            cut.push(reverse_vmapping[&ri]);
        }
    }

    cut
}


fn is_all_same(arr: &Vec<u32>) -> bool {
    arr.windows(2).all(|w| w[0] == w[1])
}


