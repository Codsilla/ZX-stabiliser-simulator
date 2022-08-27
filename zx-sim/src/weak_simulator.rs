use std::time::Duration;

use quizx::{circuit::Circuit, vec_graph::{BasisElem, GraphLike}, gate::GType::*, decompose::Decomposer};
use rand::Rng;
use crate::decompositions::Decomposition;
use quizx::vec_graph::Graph;
use crate::simulator::simulator;








pub fn qubit_by_qubit_weak_sim(c: & Circuit,decomp: fn() -> std::vec::Vec<Decomposition>) -> Vec<BasisElem> {

    let qs = c.num_qubits();
    let mut sample = Vec::with_capacity(qs);
    let pos = vec![BasisElem::Z0,BasisElem::Z1];

    for i in 0..qs {
        println!("doing qubit {}",i);
        let mut unnormalized = vec![0.0,0.0];
        for j in 0..2 {

            let mut g: Graph = c.to_graph();

            // |00..0> as input
            g.plug_inputs(&vec![BasisElem::Z0; qs]);

            g.plug_outputs(&sample);
            g.plug_outputs(&vec![pos[j].clone()]);

            // compute norm as <psi|psi>. Doubles T-count!
            g.plug(&g.to_adjoint());


            

            unnormalized[j] = simulator(&g,decomp).float_value().norm()

        }


        let mut rng = rand::thread_rng();

        if rng.gen_bool(unnormalized[0]/unnormalized.iter().sum::<f64>()) { 
            sample.push(BasisElem::Z0)
        } else {
            sample.push(BasisElem::Z1)
        }
        
    }

    sample
}


pub fn quizx_qubit_by_qubit_weak_sim(c: & Circuit,timeout:Duration) -> Vec<BasisElem> {

    let qs = c.num_qubits();
    let mut sample = Vec::with_capacity(qs);
    let pos = vec![BasisElem::Z0,BasisElem::Z1];

    for i in 0..qs {
        println!("doing qubit {}",i);
        let mut unnormalized = vec![0.0,0.0];
        for j in 0..2 {

            let mut g: Graph = c.to_graph();

            // |00..0> as input
            g.plug_inputs(&vec![BasisElem::Z0; qs]);

            g.plug_outputs(&sample);
            g.plug_outputs(&vec![pos[j].clone()]);

            // compute norm as <psi|psi>. Doubles T-count!
            g.plug(&g.to_adjoint());


            let mut d = Decomposer::new(&g)
            .with_timeout(timeout);
            d.use_cats(true);
            d.with_full_simp();
            d.decomp_all();
            d.scalar.float_value().norm();

            unnormalized[j] = d.scalar.float_value().norm();


        }


        let mut rng = rand::thread_rng();

        if rng.gen_bool(unnormalized[0]/unnormalized.iter().sum::<f64>()) { 
            sample.push(BasisElem::Z0)
        } else {
            sample.push(BasisElem::Z1)
        }
        
    }

    sample
}






pub fn gate_by_gate_weak_sim(c: & Circuit,decomp: fn() -> std::vec::Vec<Decomposition>) -> Vec<BasisElem> {

    let mut sample = vec![BasisElem::Z0; c.num_qubits()];
    let mut c_acc = Circuit::new(c.num_qubits());
    //let mut n_gates = 0;

    let mut counter = - (c.num_qubits() as isize);

    for g in &c.gates {
        //n_gates +=1;
        //println!("gates : {} / {}",n_gates,c.num_gates());
        //println!("gates : {:?}",g.t);
        c_acc.push_back(g.clone());
        

        match g.t {
            Z | CZ | CCZ | S | T | Sdg | Tdg | ZPhase => (),

            CNOT => {
                match sample[g.qs[0]] {
                    BasisElem::Z1 => sample[g.qs[1]] = sample[g.qs[1]].flipped(),
                    _ => ()
                }
            },
            //untested
            TOFF => { 
                match (sample[g.qs[0]],sample[g.qs[1]]){
                    (BasisElem::Z1,BasisElem::Z1) => sample[g.qs[2]] = sample[g.qs[2]].flipped(),
                    _ => ()
                }
            },
            //untested
            SWAP => (sample[g.qs[0]], sample[g.qs[1]]) = (sample[g.qs[1]], sample[g.qs[0]]),

            _ => {
                counter +=1 ;
                if counter % 20 == 0 {
                    println!("doing H number # {}",counter);
                }
                sample = partial_resample(&c_acc,decomp, sample, &g.qs)}
        }


    }

    sample
}


pub fn quizx_gate_by_gate_weak_sim(c: & Circuit,timeout:Duration) -> Vec<BasisElem> {

    let mut sample = vec![BasisElem::Z0; c.num_qubits()];
    let mut c_acc = Circuit::new(c.num_qubits());
    //let mut n_gates = 0;

    let mut counter = - (c.num_qubits() as isize);

    for g in &c.gates {
        //n_gates +=1;
        //println!("gates : {} / {}",n_gates,c.num_gates());
        //println!("gates : {:?}",g.t);
        c_acc.push_back(g.clone());
        

        match g.t {
            Z | CZ | CCZ | S | T | Sdg | Tdg | ZPhase => (),

            CNOT => {
                match sample[g.qs[0]] {
                    BasisElem::Z1 => sample[g.qs[1]] = sample[g.qs[1]].flipped(),
                    _ => ()
                }
            },
            //untested
            TOFF => { 
                match (sample[g.qs[0]],sample[g.qs[1]]){
                    (BasisElem::Z1,BasisElem::Z1) => sample[g.qs[2]] = sample[g.qs[2]].flipped(),
                    _ => ()
                }
            },
            //untested
            SWAP => (sample[g.qs[0]], sample[g.qs[1]]) = (sample[g.qs[1]], sample[g.qs[0]]),

            _ => {
                counter +=1 ;
                if counter % 20 == 0 {
                    println!("doing H number # {}",counter);
                }
                sample = quizx_partial_resample(&c_acc,timeout, sample, &g.qs)}
        }


    }

    sample
}





fn sample_prob_distr(distr : &[f64]) ->usize{
    

    let sum: f64 = distr.iter().sum();
    if sum < 0.999 || sum > 1.001
    {
        panic!("a probability distribution needs to sum to 1! (± 0.01 to deal with float approximation),
         but yours was {}",sum);
    }

    let mut rng = rand::thread_rng();
    let mut rando = rng.gen::<f64>();

    for (i,e) in distr.iter().enumerate() {
        if e < &rando {
            rando -= e;
        }
        else{
            return i;
        }
    }
 0
}

//implementation of https://arxiv.org/abs/2112.08499 page 2 Algortihm 2 lines 4-5
fn partial_resample(c : &Circuit,decomp: fn() -> std::vec::Vec<Decomposition>, mut sample : Vec<BasisElem>, qb_resample: &[usize]) -> Vec<BasisElem> {

    // p_t(x) for x in S (in the paper) (y_a is implied)
    let mut unormalized_distr = vec![0.0;1<<qb_resample.len()];
    for i in 0..unormalized_distr.len() {
        
        let values_sim = int_to_base_element(i,qb_resample.len());

        for (j,&q) in qb_resample.iter().enumerate() {

            sample[q] =  values_sim[j];
        }

        let mut g: Graph = c.to_graph();
        g.plug_inputs(&vec![BasisElem::Z0; c.num_qubits()]);
        g.plug_outputs(&sample);
        unormalized_distr[i] = simulator(&g,decomp).float_value().norm_sqr()
    }
 


    let renom_cst: f64 = unormalized_distr.iter().sum();
    let distr: Vec<_> = unormalized_distr.iter().map(|x| x/renom_cst).collect();
    let resampled_outpouts = int_to_base_element(sample_prob_distr(&distr),qb_resample.len());

    for (j,q) in qb_resample.iter().enumerate() {
        sample[q.clone()] =  resampled_outpouts[j];
    }

    sample
}



//implementation of https://arxiv.org/abs/2112.08499 page 2 Algortihm 2 lines 4-5
fn quizx_partial_resample(c : &Circuit,timeout:Duration, mut sample : Vec<BasisElem>, qb_resample: &[usize]) -> Vec<BasisElem> {

    // p_t(x) for x in S (in the paper) (y_a is implied)
    let mut unormalized_distr = vec![0.0;1<<qb_resample.len()];
    for i in 0..unormalized_distr.len() {
        
        let values_sim = int_to_base_element(i,qb_resample.len());

        for (j,&q) in qb_resample.iter().enumerate() {

            sample[q] =  values_sim[j];
        }

        let mut g: Graph = c.to_graph();
        g.plug_inputs(&vec![BasisElem::Z0; c.num_qubits()]);
        g.plug_outputs(&sample);

        let mut d = Decomposer::new(&g)
        .with_timeout(timeout);
        d.use_cats(true);
        d.with_full_simp();
        d.decomp_all();
        d.scalar.float_value().norm();

        unormalized_distr[i] = d.scalar.float_value().norm_sqr()
    }
 


    let renom_cst: f64 = unormalized_distr.iter().sum();
    let distr: Vec<_> = unormalized_distr.iter().map(|x| x/renom_cst).collect();
    let resampled_outpouts = int_to_base_element(sample_prob_distr(&distr),qb_resample.len());

    for (j,q) in qb_resample.iter().enumerate() {
        sample[q.clone()] =  resampled_outpouts[j];
    }

    sample
}


fn int_to_base_element(index: usize, nbqubits : usize) -> Vec<BasisElem>{

    //let elements = format!("{:b}", index);
    if index >= 1<<nbqubits
    {
        panic!("Index out of range");
    }
    format!("{:0width$b}", index, width=nbqubits).chars().map(|c| match c {
        '0' => BasisElem::Z0,
        '1' => BasisElem::Z1,
        _ => unreachable!()
    }).collect()
}


