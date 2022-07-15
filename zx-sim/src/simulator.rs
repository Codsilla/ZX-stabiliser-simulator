
use quizx::hash_graph::GraphLike;
use quizx::vec_graph::Graph;
use quizx::decompose::Decomposer;
use quizx::scalar::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use crate::decompositions::*;
use crate::utilities::to_subcomponent;


pub fn simulator(g: &Graph,decomp: fn() -> std::vec::Vec<Decomposition>) -> Scalar<Vec<isize>>{
    //Doing the first call differently saves a lot of clones
    let g = g.clone();
    simulator_internal(g,decomp,1)
}


fn simulator_internal(mut g: Graph, decomp_gen: fn() -> std::vec::Vec<Decomposition>,nb_terms:usize)-> Scalar<Vec<isize>> {
    quizx::simplify::full_simp(&mut g);

    //println!("nb term {}",nb_terms);

    if g.scalar().is_zero() {
        return ScalarN::zero();
    }

    if g.tcount() < 35 {
        let mut d = Decomposer::new(&g);
        d.use_cats(true);
        d.with_full_simp();
        d.decomp_all();
        return d.scalar;
    }

    let mut decomp = decomp_gen();

    //try all decomposition and get their alpha
    let mut inst = vec![];
    let mut alphas = vec![];

    // for mut d in decomp.clone(){
    //     let mut ng = g.clone();
    //     let verts = d.find_instance(&ng);
    //     inst.push(verts.clone());
    //     (d.to_normal_form)(&mut ng, &verts);
    //     alphas.push(d.alpha(&ng,&verts));
    // }
    for i in 0..decomp.len() {

        
        let verts = decomp[i].find_instance(&g);
        inst.push(verts.clone());

        //if no instance was found
        if decomp[i].nb_terms==usize::MAX {
            alphas.push(f64::MAX);
            continue;
        }

        if nb_terms > decomp[i].compute_alpha_until {
            alphas.push((decomp[i].approx_alpha)(&g,&verts));
        
        }else{
            let mut ng = g.clone();
            (decomp[i].to_normal_form)(&mut ng, &verts);
            alphas.push(decomp[i].alpha(&ng,&verts));
        }

    }


    let (best_decomp_index,_best_alpha) = alphas.into_iter().enumerate().reduce(|(i,x),(j,y)|{
        if x <= y { (i,x) } else { (j,y) }
    }).unwrap();

    if best_decomp_index==8 {
        println!("used Star6!!!!!")
    }

    let best_decomp = decomp[best_decomp_index];




    (best_decomp.to_normal_form)(&mut g,&inst[best_decomp_index]);

    (0..best_decomp.nb_terms)
    .into_par_iter()
    .map(|i| {
        let z = (best_decomp.get_term)(&g, &inst[best_decomp_index], &(i as usize));
        
        if z.component_vertices().len() == 0 {
            return z.scalar().clone();
        }

        let mut local_result = ScalarN::one();

        //compute scalar of subdiagrams and multiply them
        for mut subzx in to_subcomponent(&z).into_iter() {
            quizx::simplify::full_simp(&mut subzx);
            let temp = simulator_internal(subzx,decomp_gen,nb_terms*best_decomp.nb_terms);

            local_result *= temp;
            if local_result.is_zero() {
                break;
            }
        }

        //result = result + local_result;
        local_result
    })
    .reduce(|| ScalarN::zero(), |a, b| a + b)
}

