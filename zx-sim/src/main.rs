#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(unused_mut)]



mod kahypar_decomposition;
//mod kahip_cut;
mod utilities;
mod simulator;
mod decompositions;
mod weak_simulator;

use num::{Rational, Zero};
use quizx::circuit;
use quizx::decompose::Decomposer;
use quizx::hash_graph::{BasisElem, GraphLike, VType};
use quizx::scalar::ScalarN;
use quizx::vec_graph::Graph;
use rand::SeedableRng;
use rand::prelude::StdRng;
use zx_sim::utilities::hidden_shift_constructor;
use std::time::Instant;
use crate::decompositions::{Decomposition, approx_alpha_after_cut,find_star_k};
//use crate::kahip_cut::kahip_cut_finder;
use crate::kahypar_decomposition::kaHyPar_cut_finder;
use crate::simulator::simulator;
use crate::utilities::*;
use crate::weak_simulator::{qubit_by_qubit_weak_sim, gate_by_gate_weak_sim};
use std::env;

// use zx_sim::simulator::*;
// use zx_sim::utilities::*;
// use zx_sim::decompositions::Decomposition;

fn main() {

    //env::set_var("RUST_BACKTRACE", "1");
    //rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    let qs = 19;
    let seed = 3135158;

    // let c = random_iqp(qs,seed);
    // let c = hidden_shift_constructor(qs, qs, seed);
    
    
    
    let depth = 60; //nice behaviour at 130 
    let min_weight= 2 ;
    let max_weight = 4;
    
    for depth in (1..60).step_by(3){

    
        let c = random_pauli_exp(qs, depth, seed, min_weight, max_weight);

        let mut g: Graph = c.clone().to_graph();

        g.plug_inputs(&vec![BasisElem::Z0; qs]);
        g.plug_outputs(&vec![BasisElem::Z0; qs]);
        

        
        println!("T-count before simplification {}",g.tcount());
        quizx::simplify::full_simp(&mut g);
        println!("T-count after simplification {}",g.tcount());
        println!("_________________");

        let time = Instant::now();
        println!("got       {}",simulator(&g,Decomposition::without_cuts).to_float());
        println!("Simulator time without_cuts  : {:.2?}",time.elapsed());

        let time = Instant::now();
        let mut d = Decomposer::new(&g);
        d.use_cats(true);
        d.with_full_simp();
        d.decomp_all();
        println!("should be {}", d.scalar.to_float());
        println!("Quizx Simulator time  : {:.2?}",time.elapsed());

    }

    //let mut decomp = Decomposition::all_decomp();


    

    // let mut inst = vec![];
    // let mut alphas = vec![];
    // for i in 0..decomp.len() {

    //     let mut ng = g.clone();
    //     let verts = decomp[i].find_instance(&ng);
    //     inst.push(verts.clone());
    //     (decomp[i].to_normal_form)(&mut ng, &verts);
    //     alphas.push(decomp[i].alpha(&ng,&verts));
    // }
    // dbg!(alphas);

    //let mut decomp = Decomposition::without_cuts();
    // let inst = decomp.find_instance(&g);

    // dbg!(g.component_vertices().len());
    // //println!("{}",g.to_dot());
    // println!("-------------------");
    // println!("cut size {}",inst.len());
    // println!("The cut {:?}",inst);
    // println!("nb term {}",decomp.nb_terms.clone());

    // let mut comp = (decomp.get_term)(&g,&inst,&0);
    // dbg!(comp.component_vertices().len());
    
    // //println!("{}",comp.to_dot());
    
    // quizx::simplify::full_simp(&mut comp);
    // dbg!(comp.component_vertices().len());
    // println!("T count after : {}",comp.tcount());
    // //println!("{}",comp.to_dot());
    
    // let mut lastcomp = (decomp.get_term)(&g,&inst,&(1<<inst.len()-1));
    // quizx::simplify::full_simp(&mut lastcomp);
    // dbg!(lastcomp.component_vertices().len());
    // println!("T count after  last : {}",lastcomp.tcount());

    // println!("--------------Kahypar cut------------------------");
    // println!("cut size : {}",inst[2].len());
    // //println!("The cut {:?}",inst[2]);
    // println!("alpha {:?}",approx_alpha_after_cut(&g,&inst[2]));

    // println!("--------------unbalanced Kahypar cut------------------------");
    // let cut = kaHyPar_cut_finder(&g,0.9);
    // println!("cut size : {}",cut.len());
    // //println!("The cut {:?}",inst[2]);
    // println!("alpha {:?}",approx_alpha_after_cut(&g,&cut));


    // println!("--------------Kahip cut------------------------");
    // let kahipcut = kahip_cut_finder(&g,0.1);
    // println!("cut size : {}",kahipcut.len());
    // //println!("The cut {:?}",kahipcut);
    // println!("alpha {:?}",approx_alpha_after_cut(&g,&kahipcut));

    // println!("--------------unbalanced cut------------------------");
    // let kahipcut = kahip_cut_finder(&g,0.9);
    // println!("cut size : {}",kahipcut.len());
    // //println!("The cut {:?}",kahipcut);
    // println!("alpha {:?}",approx_alpha_after_cut(&g,&kahipcut));


    // let time = Instant::now();
    // println!("got       {}",simulator(&g,Decomposition::without_cuts).to_float());
    // println!("Simulator time without_cuts  : {:.2?}",time.elapsed());

    // let time = Instant::now();
    // let mut d = Decomposer::new(&g);
    // d.use_cats(true);
    // d.with_full_simp();
    // d.decomp_all();
    // println!("should be {}", d.scalar.to_float());
    // println!("Quizx Simulator time  : {:.2?}",time.elapsed());


    // let time = Instant::now();
    // println!("got       {}",simulator(&g,Decomposition::with_best_t).to_float());
    // println!("Simulator time with_best_t  : {:.2?}",time.elapsed());

    // println!("Strong sim");

    // let time = Instant::now();
    // println!("got       {}",simulator(&g,Decomposition::without_cuts).to_float());
    // println!("Simulator time all_decomp  : {:.2?}",time.elapsed());

    // println!("Weak sim");

    // let time = Instant::now();
    // println!(" {:?}" ,gate_by_gate_weak_sim(&c, Decomposition::without_cuts));
    // println!("took  : {:.2?}",time.elapsed());


    // let time = Instant::now();
    // println!(" {:?}" ,qubit_by_qubit_weak_sim(&c, Decomposition::all_decomp));
    // println!("took  : {:.2?}",time.elapsed());




}






//Test cats and stars

// println!("----------cats-----------");

// for i in 7..=10{
    //     let mut g = make_cat_n(i,false);
    
    //     //dbg!(cat_distribution(&g));
//     let time = Instant::now();
//     println!("got       {}",simulator(&g,Decomposition::test_plz).to_float());
//    // println!("Simulator time without_cuts  : {:.2?}",time.elapsed());

    
//     let time = Instant::now();
//     let mut d = Decomposer::new(&g);
//     d.use_cats(true);
//     d.with_full_simp();
//     d.decomp_all();
//     println!("should be {}", d.scalar.to_float());
//     //println!("Quizx Simulator time  : {:.2?}",time.elapsed());

// }

// println!("---------Stars-----------");

// for i in 2..=10{

//     if i !=2 && i!=6 && i!=10 {
//         continue;
//     }
//     let mut g = make_star_n(i,false);

//     //dbg!(cat_distribution(&g));
//     let time = Instant::now();
//     println!("got       {}",simulator(&g,Decomposition::test_plz).to_float());
//    // println!("Simulator time without_cuts  : {:.2?}",time.elapsed());

    
//     let time = Instant::now();
//     let mut d = Decomposer::new(&g);
//     d.use_cats(true);
//     d.with_full_simp();
//     d.decomp_all();
//     println!("should be {}", d.scalar.to_float());
//     //println!("Quizx Simulator time  : {:.2?}",time.elapsed());

// }