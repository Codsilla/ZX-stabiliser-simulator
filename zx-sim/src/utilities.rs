use num::{One, Rational};
use quizx::circuit::Circuit;
use quizx::vec_graph::Graph;
use rand::Rng;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::hash::Hash;
use quizx::hash_graph::{GraphLike, VType};


pub fn cat_distribution(g: &Graph) -> Vec<usize>{

    let mut cats:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer()).map(|x| g.neighbors(x).len()).collect();

    cats.sort();

    cats
}

//computes the value of alpha given the resulting t-count 
pub fn alpha_from_after(g: &Graph, tcount_after: isize, nb_terms: usize) -> f64 {

    let t_removed = g.tcount() as isize - tcount_after;
    (nb_terms as f64).log2() / (t_removed as f64)
}


pub fn to_subcomponent(g :& Graph) -> Vec<Graph>{


    let comp = g.component_vertices();
    let mut subzx = vec![g.clone();comp.len()];


    for i in 0..comp.len(){
        for v in g.vertices() {

            if !comp[i].contains(&v) {
                subzx[i].remove_vertex(v);
            }
        }

        //remove phase except for the first one
        if i!=0 {
            subzx[i].scalar_mut().set_one();

        }
    }

    subzx
}






// from https://stackoverflow.com/questions/64262297/rust-how-to-find-n-th-most-frequent-element-in-a-collection
pub fn most_frequent<T>(array: &[T], k: usize) -> (usize, &T) where T: Hash + Eq + Ord,{
    let mut map = HashMap::with_capacity(array.len());
    for x in array {
        *map.entry(x).or_default() += 1;
    }

    let mut heap = BinaryHeap::with_capacity(k + 1);
    for (x, count) in map.into_iter() {
        if heap.len() < k {
            heap.push(Reverse((count, x)));
        } else {
            let &Reverse((min, _)) = heap.peek().unwrap();
            if min < count {
                heap.pop();
                heap.push(Reverse((count, x)));
            }
        }
    }
    let result:Vec<(usize, &T)> = heap.into_sorted_vec().into_iter().map(|r| r.0).collect();
    result[k - 1]
}


pub fn make_star_n(n: usize)->Graph{
    let mut g = Graph::new();
    g.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
    
    for i in 1..=n{
        g.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
        g.add_edge_smart(0,i,quizx::hash_graph::EType::H);
    }
    g
}



pub fn hidden_shift_constructor(qs :usize, n_ccz : usize, seed: usize) ->Circuit {


    let debug = true;
    if debug { println!("qubits: {}, # ccz: {}, seed: {}", qs, n_ccz, seed); }

    
    let (c,_shift) = Circuit::random_hidden_shift()
        .qubits(qs)
        .n_ccz(n_ccz) // T = CCZ * 2 * 7
        .seed((seed*qs*n_ccz) as u64)
        .build();

    c
}


pub fn random_iqp(nqubits: usize) -> Circuit {


    let mut c_acc = Circuit::new(nqubits);
    //fisrt layer of Hadamard + random T
    for i in 0..nqubits{
        c_acc.add_gate( "h", vec![i]);
        let t_phase = rand::thread_rng().gen_range(0..8);
        c_acc.add_gate_with_phase("rz", vec![i], Rational::new(t_phase,4));
    }

    for i in 0..nqubits{
        for j in 0..nqubits{
            if j==i {continue} // the parser did not crash with cx [i,i]!!!!!

            let n_sgate = rand::thread_rng().gen_range(0..3);
            
            for _ in 0..n_sgate {
                //implementation of a CS gate
                c_acc.add_gate_with_phase("rz", vec![i], Rational::new(1,4));
                c_acc.add_gate_with_phase("rz", vec![j], Rational::new(1,4));
                c_acc.add_gate("cx", vec![i,j]);
                c_acc.add_gate_with_phase("rz", vec![j], Rational::new(7,4));
                c_acc.add_gate("cx", vec![i,j]);
            }
        }
        //final layer of Hadamard
        c_acc.add_gate( "h", vec![i]);
    }



    c_acc
}