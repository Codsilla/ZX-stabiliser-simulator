use num::{One, Rational};
use quizx::circuit::Circuit;
use quizx::scalar::{ScalarN, Sqrt2};
use quizx::vec_graph::Graph;
use rand::{Rng, SeedableRng};
use rand::prelude::StdRng;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::hash::Hash;
use quizx::hash_graph::{GraphLike, VType};
use std::collections::HashSet;


pub fn cat_distribution(g: &Graph) -> Vec<usize>{

    let mut cats:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer()).map(|x| g.neighbors(x).len()).collect();

    cats.sort();

    cats
}


pub fn star_distribution(g: &Graph) -> Vec<usize>{

    let potential_stars:Vec<usize> = g.vertices().filter(|&v| *g.phase(v).denom() == 4).collect();

    let mut stars_degree = vec![];

    'outer: for v in potential_stars {
        for neigh in g.neighbors(v) {
            if *g.phase(neigh).denom() != 4 { continue 'outer;}
        }

        stars_degree.push(g.neighbors(v).len()); 
    }

    stars_degree.sort();

    stars_degree
}

pub fn find_stars_center(g: &Graph) -> Vec<usize>{

    let potential_stars:Vec<usize> = g.vertices().filter(|&v| *g.phase(v).denom() == 4).collect();

    let mut stars = vec![];

    'outer: for v in potential_stars {
        for neigh in g.neighbors(v) {
            if *g.phase(neigh).denom() != 4 { continue 'outer;}
        }

        stars.push(v); 
    }

    stars

}


pub fn star_differ_distribution(g: &Graph) -> Vec<usize>{

    let stars = find_stars_center(g);
    let mut count = vec![0; 10];

    for i in 0..stars.len() {

        let ni:HashSet<usize> = HashSet::from_iter(g.neighbors(stars[i]));

        for j in i+1..stars.len() {

            if ni.contains(&stars[j]) { continue  }

            let nj:HashSet<usize> = HashSet::from_iter(g.neighbors(stars[j]));
            
            let symdif = ni.symmetric_difference(&nj).count();

            if symdif < 10 {
                count[symdif] += 1;
            }

        }
    }

    count
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


pub fn make_star_n(n: usize, add_pi_center: bool)->Graph{
    let mut g = Graph::new();
    g.add_vertex_with_phase(VType::Z,Rational::new(1, 4));

    if add_pi_center {g.add_to_phase(0, Rational::new(1, 1))}
    
    for i in 1..=n{
        g.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
        g.add_edge_smart(0,i,quizx::hash_graph::EType::H);
    }
    g
}

pub fn make_cat_n(n: usize, add_pi_center: bool)->Graph{
    let mut g = Graph::new();
    g.add_vertex_with_phase(VType::Z,Rational::new(0, 1));

    if add_pi_center {g.add_to_phase(0, Rational::new(1, 1))}
    
    for i in 1..=n{
        g.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
        g.add_edge_smart(0,i,quizx::hash_graph::EType::H);
    }
    g
}




pub fn hidden_shift_constructor(qs :usize, n_ccz : usize, seed: usize) ->Circuit {



    println!("qubits: {}, # ccz: {}, seed: {}", qs, n_ccz, seed);

    
    let (c,_shift) = Circuit::random_hidden_shift()
        .qubits(qs)
        .n_ccz(n_ccz) // T = CCZ * 2 * 7
        .seed((seed*qs*n_ccz) as u64)
        .build();

    c
}


pub fn random_iqp(nqubits: usize,seed : u64) -> Circuit {

    let mut rng =StdRng::seed_from_u64(seed);
    let mut c_acc = Circuit::new(nqubits);
    //fisrt layer of Hadamard + random T
    for i in 0..nqubits{
        c_acc.add_gate( "h", vec![i]);
        let t_phase = rng.gen_range(0..8);
        c_acc.add_gate_with_phase("rz", vec![i], Rational::new(t_phase,4));
    }

    for i in 0..nqubits{
        for j in i+1..nqubits{
            if j==i {continue} // the parser did not crash with cx [i,i]!!!!!

            let n_sgate = rng.gen_range(0..=3);

            
            if n_sgate == 0 {continue;}
        //implementation of a power of a CS gate
           c_acc.add_gate_with_phase("rz", vec![i], Rational::new(n_sgate,4));
           c_acc.add_gate_with_phase("rz", vec![j], Rational::new(n_sgate,4));
           c_acc.add_gate("cx", vec![i,j]);
           c_acc.add_gate_with_phase("rz", vec![j], Rational::new(-n_sgate,4));
           c_acc.add_gate("cx", vec![i,j]);
        
        }
        //final layer of Hadamard
        c_acc.add_gate( "h", vec![i]);
    }


    //for i in 0..nqubits{
    //    c_acc.add_gate( "h", vec![i]);
    //}


    c_acc
}

// let n_sgate = rng.gen_range(0..4);
// if n_sgate == 0 {continue;}
// // for _ in 0..n_sgate {
//     //implementation of a CS gate
//     c_acc.add_gate_with_phase("rz", vec![i], Rational::new(n_sgate,4));
//     c_acc.add_gate_with_phase("rz", vec![j], Rational::new(n_sgate,4));
//     c_acc.add_gate("cx", vec![i,j]);
//     c_acc.add_gate_with_phase("rz", vec![j], Rational::new(-n_sgate,4));
//     c_acc.add_gate("cx", vec![i,j]);
// //}


// if n_sgate == 0 {continue;}
// for _ in 0..=n_sgate {
//     //implementation of a CS gate
//     c_acc.add_gate_with_phase("rz", vec![i], Rational::new(1,4));
//     c_acc.add_gate_with_phase("rz", vec![j], Rational::new(1,4));
//     c_acc.add_gate("cx", vec![i,j]);
//     c_acc.add_gate_with_phase("rz", vec![j], Rational::new(-1,4));
//     c_acc.add_gate("cx", vec![i,j]);
// }


pub fn random_pauli_exp(qs :usize, depth : usize,seed: u64, min_weight: usize, max_weight: usize) ->Circuit  {
    
    println!("qubits: {}, depth: {}, min_weight: {}, max_weight: {}, seed: {}",qs, depth, min_weight, max_weight, seed);
    
    let c = Circuit::random_pauli_gadget()
    .qubits(qs)
    .depth(depth)
    .seed(seed)
    .min_weight(min_weight)
    .max_weight(max_weight)
    .build();
    

    c
}

//buggy sometimes give the error 'Parallel edges only supported between Z and X vertices'
//example fix_t_expectation_random_circuit(40,400,30,1355); after simplification
pub fn fix_t_expectation_random_circuit(qs:usize,depth : usize, nb_t : usize, seed:u64)->Circuit {


    let c = Circuit::random().seed(seed)
                .qubits(qs)
                .depth(depth)
                .p_t((nb_t as f32)/(depth as f32))
                .with_cliffords()
                .build();

    c
}



//There is a bug here apparently (this function is not used anywhere else)
//fuse out a t and cut its wire (for catification) 
pub fn remove_t_cut_wire(g: & Graph,spider : usize, index: usize) -> Graph{
    
    if index > 1 {panic!("Look for the index {} of the remove_t_cut_wire",index)}

    let mut gn = g.clone();
    let reso = gn.add_vertex_with_phase(VType::Z,Rational::new(index as isize, 1));
    let fuseout = gn.add_vertex_with_phase(VType::Z,gn.phase(spider));
    gn.add_edge(reso, fuseout);
    gn.set_phase(spider, Rational::new(index as isize,1));

    //normalization
    *gn.scalar_mut() *= ScalarN::sqrt2_pow(-2);
    gn
}



//fuse out a cat into two by adding T gates in the centere (for cat splits)
pub fn cat_spliter(g: & Graph, center : usize, k : usize, n : usize) -> (Graph,Vec<usize>,Vec<usize>){
    
    if g.degree(center) != k+n-2 {panic!("Tried to split a cat {} into a cat {} and a cat {}",g.degree(center),k,n)}


    let mut gn = g.clone();

    let catn_center = gn.add_vertex(VType::Z);

    let mut catk = vec![center];
    let mut catn = vec![catn_center];

    let neigh = gn.neighbor_vec(center);

    //unfuse
    for i in 0..(n-1) {
        
        gn.add_edge_with_type(neigh[i], catn_center, quizx::hash_graph::EType::H);
        gn.remove_edge(neigh[i], center);

    }

    //add new Ts
    let newt = gn.add_vertex_with_phase(VType::Z,Rational::new(1,4));
    let new_anti_t = gn.add_vertex_with_phase(VType::Z,Rational::new(-1,4));
    gn.add_to_phase(new_anti_t, gn.phase(center));
    gn.set_phase(center, Rational::new(0,1));

    gn.add_edge_with_type(newt,new_anti_t,quizx::hash_graph::EType::N);

    gn.add_edge_with_type(newt,catn[0],quizx::hash_graph::EType::H);
    gn.add_edge_with_type(new_anti_t,catk[0],quizx::hash_graph::EType::H);

    catn.append(&mut gn.neighbor_vec(catn[0]));
    catk.append(&mut gn.neighbor_vec(catk[0]));



    (gn,catk,catn)
}