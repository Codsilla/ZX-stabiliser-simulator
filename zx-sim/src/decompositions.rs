

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(non_snake_case)]


use num::*;
use array_tool::vec::Intersect;
use quizx::hash_graph::{BasisElem, GraphLike, VType, EType};
use quizx::scalar::*;
use quizx::vec_graph::{self, Graph};
use std::cmp::{Eq, Ord, Reverse};
use std::cmp;


use crate::kahypar_decomposition::kaHyPar_cut_finder;
use crate::utilities::*;

const THRESHOLD_COMPUTE_ALPHA: usize = 30;


#[derive(Clone, Copy)]
pub struct Decomposition {
    finder: fn(&Graph) -> Vec<usize>, 
    pub nb_terms: usize, //if not known at compile time (e.g. for cuts) put 0
    pub to_normal_form: fn(&mut Graph,&Vec<usize>),
    pub get_term: fn(&Graph,&Vec<usize>,&usize) -> Graph,
    pub approx_alpha: fn(&Graph,&Vec<usize>) -> f64,
    pub compute_alpha_until: usize
}


impl Decomposition {

    pub fn find_instance(&mut self,g:&Graph)-> Vec<usize>{
        let inst = (self.finder)(g);
        
        //for cuts
        if self.nb_terms == 0 {
            self.nb_terms = 1 << inst.len();
        }

        if inst.len() == 0 {
            self.nb_terms = usize::MAX;
        }

        inst
    }



    //give the score of a decomposition (effective alpha)
    pub fn alpha(&self,g: &Graph,verts: &Vec<usize>) -> f64 {

        if self.nb_terms == usize::MAX {
            return f64::MAX;
        }

        // if g.tcount() < self.compute_alpha_until {
        //     return (self.approx_alpha)(&g,&verts);
        // }

        let mut g = g.clone();
        (self.to_normal_form)(&mut g,verts);


        let nb_term_check = if self.nb_terms>5 {1} else{self.nb_terms};//

        let mut tcounts : Vec<isize> = Vec::new();
        
        let base_alpha = (self.approx_alpha)(&g,&verts);
        //check max T in connected component for the last few terms (last because better heuristic for cuts)
        for i in self.nb_terms-nb_term_check..self.nb_terms {
            
            let mut d = (self.get_term)(&g,&verts,&i);   
            quizx::simplify::full_simp(&mut d);

            let comp = d.component_vertices();

            if comp.len()==0 {
                tcounts.push(0);
                continue;
            }
            
            //count the T in each connected component
            for mut c in comp {
                c.sort();
                c.dedup();
                let t = d.tcount();
                
        
                for i in 0..c.len() {
                    let v = c[i];
                    d.remove_vertex(v);
                }
                let t_component = t as isize - d.tcount() as isize;
                tcounts.push(t_component);

                //helps to compute faster
                //let alpha_comp = (self.nb_terms as f64).log2() / ((g.tcount() as isize - t_component) as f64);
                if alpha_from_after(&g,t_component,self.nb_terms) + 0.1 > base_alpha {
                    return base_alpha;
                }
            }
            
        }
        //println!("t count term {:?}", tcounts);
        let max_t_after = *tcounts.iter().max().expect("something went wrong");

        let t_removed = g.tcount() - max_t_after as usize;


        //dbg!(t_removed);
        //dbg!((self.nb_terms as f64).log2());
        //let alpha = (self.nb_terms as f64).log2() / (t_removed as f64);
        let alpha = alpha_from_after(&g,max_t_after,self.nb_terms);
        // if self.nb_terms > 10 {
        //     // if alpha < 0.24{
        //     //     return 1.0;
        //     // }
        //     let alpha = alpha + 0.00;
        // }
        return alpha
    }



    pub fn all_decomp()->Vec<Decomposition>{
        vec![Decomposition::trivial_decomp(),Decomposition::trivial_in_cat3_decomp(),Decomposition::cut_decomp()
        ,Decomposition::magic5_2_decomp(),Decomposition::cat3_decomp(),Decomposition::cat4_decomp()
        ,Decomposition::cat5_decomp(),Decomposition::cat6_decomp()]
    }

    pub fn without_cuts()->Vec<Decomposition>{
        vec![Decomposition::trivial_decomp(),Decomposition::trivial_in_cat3_decomp()
        ,Decomposition::magic5_2_decomp(),Decomposition::cat3_decomp(),Decomposition::cat4_decomp()
        ,Decomposition::cat5_decomp(),Decomposition::cat6_decomp()]
    }


    pub fn without_trivial_and_cuts()->Vec<Decomposition>{
        vec![Decomposition::trivial_in_cat3_decomp()
        ,Decomposition::magic5_2_decomp(),Decomposition::cat3_decomp(),Decomposition::cat4_decomp()
        ,Decomposition::cat5_decomp(),Decomposition::cat6_decomp()]
    }

    pub fn with_stars()->Vec<Decomposition>{
        vec![Decomposition::trivial_decomp(),Decomposition::trivial_in_cat3_decomp(),Decomposition::cut_decomp()
        ,Decomposition::magic5_2_decomp(),Decomposition::cat3_decomp(),Decomposition::cat4_decomp()
        ,Decomposition::cat5_decomp(),Decomposition::cat6_decomp(),Decomposition::star6_decomp(),Decomposition::star3_decomp()]
    }

}



fn first_ts(g: &Graph,k : usize) -> Vec<usize> {
    let mut t = vec![];

    for v in g.vertices() {
        if *g.phase(v).denom() == 4 { t.push(v); }
        if t.len() == k { return t; }
    }

    vec![]
}


pub fn most_connected_t_incat3(g :&Graph) -> Vec<usize>{

    let cats3:Vec<usize> = g.vertices().filter(|&v| g.phase(v).is_integer() && g.degree(v) == 3).collect();
    if cats3.len() == 0 {return vec![]}
    let mut ts = vec![];

    for c in cats3 {
        ts.append(&mut g.neighbor_vec(c));
    }

    let (count,x) = most_frequent(&ts, 1);
    
    //println!("best t gate is in {} cat 3",count);
    vec![*x]
}

fn trivial_normal_form(g:&mut Graph,verts: &Vec<usize>){}

fn replace_t0(g:  &Graph, verts: &[usize]) -> Graph {
    // println!("replace_t0");
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);
    let w = g.add_vertex(VType::Z);
    g.add_edge_with_type(verts[0], w, EType::H);
    g.add_to_phase(verts[0], Rational::new(-1,4));
    g
}

fn replace_t1(g: &Graph, verts: &[usize]) -> Graph {
    // println!("replace_t1");
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1,0,1,0]);
    let w = g.add_vertex_with_phase(VType::Z, Rational::one());
    g.add_edge_with_type(verts[0], w, EType::H);
    g.add_to_phase(verts[0], Rational::new(-1,4));
    g
}

fn trivial_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize)-> Graph{

    match i {
        0 => replace_t0(g, verts),
        1 => replace_t1(g, verts),
        other => panic!("Tried to access the index {} of the trivial decomposition",other)
    }

}

//TODO: add quicker approx for trivial reaplace (count cat3-arity)

impl Decomposition {
    pub fn trivial_decomp()->Decomposition{
        Decomposition { finder: |g: &Graph| -> Vec<usize>{ first_ts(g, 2) }, nb_terms: 2, to_normal_form: trivial_normal_form, get_term: trivial_replace_index, approx_alpha: |g ,vert| 1.0, compute_alpha_until: 1  }
    }
    pub fn trivial_in_cat3_decomp()->Decomposition{
        Decomposition { finder: |g: &Graph| -> Vec<usize>{ most_connected_t_incat3(g) }, nb_terms: 2, to_normal_form: trivial_normal_form, get_term: cut_index, approx_alpha: approx_alpha_after_cut, compute_alpha_until:usize::MAX  }
    }
}

fn replace_magic5_0(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(1, vec![1, 0, 0, 0]);    
    for &v in verts {
        g.add_to_phase(v, Rational::new(-1,4));
        g.add_edge_smart(v, verts[0], EType::N);
    }
    g.add_to_phase(verts[0], Rational::new(-3,4));
    g
}

fn replace_magic5_1(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();  
    *g.scalar_mut() *= ScalarN::Exact(1, vec![-1, 0, 1, 0]);
    let p = g.add_vertex(VType::Z);
    for &v in verts {
        g.add_to_phase(v, Rational::new(-1,4));
        g.add_edge_with_type(v, p, EType::H);
    }
    let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1,4));
    g.add_edge_with_type(w, p, EType::H);
    g
}

fn replace_magic5_2(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(9, vec![0, -1, 0, 0]);
    let p = g.add_vertex(VType::Z);
    let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1,4));
    g.add_edge_with_type(p, w, EType::H);
    for i in 0..verts.len() {
        g.add_to_phase(verts[i], Rational::new(-1,4));
        g.add_edge_with_type(verts[i], p, EType::H);
        g.add_edge_with_type(verts[i], w, EType::H);
        for j in i+1..verts.len() {
            g.add_edge_smart(verts[i], verts[j], EType::H);
        }
    }
    g
}

fn magic5_2_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize) -> Graph{
    match i {
        0 => replace_magic5_0(g, verts),
        1 => replace_magic5_1(g, verts),
        2 => replace_magic5_2(g, verts),
        other => panic!("Tried to access the index {} of the magic5_2 decomposition",other)
    }
}

impl Decomposition {
   pub fn magic5_2_decomp()->Decomposition{
        Decomposition { finder: |g: &Graph| -> Vec<usize>{ first_ts(g, 5) }, nb_terms: 3, to_normal_form: trivial_normal_form, get_term: magic5_2_replace_index, approx_alpha: |g ,vert| 0.396, compute_alpha_until:THRESHOLD_COMPUTE_ALPHA  }
    }
}


//Cats

fn find_cat_k(g: &Graph,k : usize) -> Vec<usize>{

    let center = g.vertices().find(|&v| g.phase(v).is_integer() && g.degree(v) == k).unwrap_or(usize::MAX);
    
    if center == usize::MAX { 
        return Vec::new()
    }

    let mut cat = vec![center]; 
    cat.append(&mut g.neighbor_vec(center));
    cat
}

fn cat_normal_form(g : &mut Graph, verts: &Vec<usize>){
    
    if verts.len() == 0 {return;}

    if g.phase(verts[0]).numer() == &1 {
        g.set_phase(verts[0], Rational::new(0,1));
        let mut neigh = g.neighbor_vec(verts[1]);
        neigh.retain(|&x| x != verts[0]);
        for &v in &neigh{
            g.add_to_phase(v, Rational::new(1,1));
        }
        let tmp = g.phase(verts[1]);
        *g.scalar_mut() *= ScalarN::from_phase(tmp);
        g.set_phase(verts[1], g.phase(verts[1])*Rational::new(-1,1));
    }
}

fn replace_cat3_0(g: &Graph, verts: &[usize])-> Graph {
    let mut g = g.clone();
    let mut verts = Vec::from(verts);

    let w = g.add_vertex(VType::Z);
    let v = g.add_vertex(VType::Z);
    g.add_edge_with_type(v, w, EType::H);
    g.add_edge_with_type(v, verts[0], EType::H);
    verts.push(v);   

    *g.scalar_mut() *= ScalarN::Exact(0, vec![0, 0, 1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
    }
    g
}

fn replace_cat3_1(g: &Graph, verts: &[usize])-> Graph {
 
    // same as replace_cat6_0, only with a different scalar
    let mut g = g.clone();  
    let mut verts = Vec::from(verts);

    let w = g.add_vertex(VType::Z);
    let v = g.add_vertex(VType::Z);
    g.add_edge_with_type(v, w, EType::H);
    g.add_edge_with_type(v, verts[0], EType::H);
    verts.push(v);  

    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, -1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational::new(-1,2));
    g
    
}

fn cat3_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize) -> Graph{

    match i {
        0 => replace_cat3_0(g, verts),
        1 => replace_cat3_1(g, verts),
        other => panic!("Tried to access the index {} of the replace_cat4 decomposition",other)
    }

}

impl Decomposition {
   pub fn cat3_decomp()->Decomposition{
        Decomposition { finder: |g: &Graph| -> Vec<usize>{ find_cat_k(g, 3) }, nb_terms: 2, to_normal_form: cat_normal_form, get_term: cat3_replace_index, approx_alpha: |g,vert| 0.333, compute_alpha_until: THRESHOLD_COMPUTE_ALPHA  }
    }
}

fn replace_cat4_0(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(0, vec![0, 0, 1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
    }
    g
}

fn replace_cat4_1(g: &Graph, verts: &[usize]) -> Graph {
    // same as replace_cat6_0, only with a different scalar
    let mut g = g.clone();  
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, -1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational::new(-1,2));
    g
}

fn cat4_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize) -> Graph{

    match i {
        0 => replace_cat4_0(g, verts),
        1 => replace_cat4_1(g, verts),
        other => panic!("Tried to access the index {} of the replace_cat4 decomposition",other)
    }

}

impl Decomposition {
   pub fn cat4_decomp()->Decomposition{
        Decomposition { finder: |g: &Graph| -> Vec<usize>{ find_cat_k(g, 4) }, nb_terms: 2, to_normal_form: cat_normal_form, get_term: cat4_replace_index,approx_alpha: |g,vert| 0.25, compute_alpha_until: THRESHOLD_COMPUTE_ALPHA  }
    }
}

fn replace_cat5_0(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();

    let mut verts = Vec::from(verts);

    let w = g.add_vertex(VType::Z);
    let v = g.add_vertex(VType::Z);
    g.add_edge_with_type(v, w, EType::H);
    g.add_edge_with_type(v, verts[0], EType::H);
    verts.push(v); 

    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, 0, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational::new(-1,2));
    g
}

fn replace_cat5_1(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();  
    let mut verts = Vec::from(verts);

    let w = g.add_vertex(VType::Z);
    let v = g.add_vertex(VType::Z);
    g.add_edge_with_type(v, w, EType::H);
    g.add_edge_with_type(v, verts[0], EType::H);
    verts.push(v); 

    *g.scalar_mut() *= ScalarN::Exact(-1, vec![-1, 0, 1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
    }
    g
}

fn replace_cat5_2(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    let mut verts = Vec::from(verts);

    let w = g.add_vertex(VType::Z);
    let v = g.add_vertex(VType::Z);
    g.add_edge_with_type(v, w, EType::H);
    g.add_edge_with_type(v, verts[0], EType::H);
    verts.push(v); 
    
    *g.scalar_mut() *= ScalarN::Exact(7, vec![0, -1, 0, 0]);
    for i in 1..verts.len() {
        g.add_to_phase(verts[i], Rational::new(-1,4));
        for j in i+1..verts.len() {
            g.add_edge_smart(verts[i], verts[j], EType::H);
        }
    }
    g
}

fn cat5_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize) -> Graph{

    match i {
        0 => replace_cat5_0(g, verts),
        1 => replace_cat5_1(g, verts),
        2 => replace_cat5_2(g, verts),
        other => panic!("Tried to access the index {} of the replace_cat5 decomposition",other)
    }

}

impl Decomposition {
   pub fn cat5_decomp()->Decomposition{
        Decomposition { finder: |g: &Graph| -> Vec<usize>{ find_cat_k(g, 5) }, nb_terms: 3, to_normal_form: cat_normal_form, get_term: cat5_replace_index, approx_alpha: |g,vert| 0.317, compute_alpha_until: THRESHOLD_COMPUTE_ALPHA }
    }
}

fn replace_cat6_0(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, 0, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational::new(-1,2));
    g
}

fn replace_cat6_1(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();  
    *g.scalar_mut() *= ScalarN::Exact(-1, vec![-1, 0, 1, 0]);    
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational::new(-1,4));
    }
    g
}

fn replace_cat6_2(g: &Graph, verts: &[usize]) -> Graph {
    let mut g = g.clone();
    *g.scalar_mut() *= ScalarN::Exact(7, vec![0, -1, 0, 0]);
    for i in 1..verts.len() {
        g.add_to_phase(verts[i], Rational::new(-1,4));
        for j in i+1..verts.len() {
            g.add_edge_smart(verts[i], verts[j], EType::H);
        }
    }
    g
}

fn cat6_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize) -> Graph{

    match i {
        0 => replace_cat6_0(g, verts),
        1 => replace_cat6_1(g, verts),
        2 => replace_cat6_2(g, verts),
        other => panic!("Tried to access the index {} of the replace_cat6 decomposition",other)
    }

}

impl Decomposition {
   pub fn cat6_decomp()->Decomposition{
        Decomposition { finder: |g: &Graph| -> Vec<usize>{ find_cat_k(g, 6) }, nb_terms: 3, to_normal_form: cat_normal_form, get_term: cat6_replace_index, approx_alpha: |g,vert| 0.264, compute_alpha_until: THRESHOLD_COMPUTE_ALPHA }
    }
}


//Cuts

//              (g: &Graph,verts: &Vec<usize>, i: &usize)-> Graph{
fn cut_index(g :& Graph,cuts :& Vec<usize>, index : &usize) -> Graph {

    let mut ng = g.clone();
    let mut i = 0;

    for c in cuts {

        //normalization for basis state
        *ng.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);

        let nv = ng.add_vertex(VType::X);
        ng.add_edge(*c, nv);


        if index >> i & 1 == 1 {
            ng.set_phase(nv, Rational::new(1,1));
        }
        i += 1;
        
    }
    //quizx::simplify::full_simp(&mut ng);
    //println!("T count :{}",ng.tcount());
    ng
}


fn approx_t_after_cut(g :& Graph,cuts :& Vec<usize>) -> usize{

    let mut  g = g.clone();
    
    for &v in cuts {
        g.remove_vertex(v);
    }

    quizx::simplify::full_simp(&mut g);

    let comp = g.component_vertices();

    if comp.len() == 0 { return 0 }

    let mut tcounts : Vec<isize> = Vec::new();
    //count number of T in each component by removing them one by one
    for mut c in comp {
        c.sort();
        c.dedup();
        let t = g.tcount();
        

        for i in 0..c.len() {
            let v = c[i];
            g.remove_vertex(v);
        }

        tcounts.push(t as isize - g.tcount() as isize);
    }
    let t_after = *tcounts.iter().max().expect("something went wrong");

    t_after as usize
}

fn approx_alpha_after_cut(g :& Graph,cuts :& Vec<usize>) -> f64 {
    alpha_from_after(&g,approx_t_after_cut(&g,&cuts) as isize, 1<< cuts.len())
}


impl Decomposition {
    pub fn cut_decomp()->Decomposition{
         Decomposition { finder: |g: &Graph| -> Vec<usize>{ kaHyPar_cut_finder(g,0.3) }, nb_terms: 0, to_normal_form: trivial_normal_form, get_term: cut_index, approx_alpha: approx_alpha_after_cut, compute_alpha_until:THRESHOLD_COMPUTE_ALPHA }
     }
 }



 //Stars

 fn find_star_k(g: &Graph,k : usize) -> Vec<usize>{

    let t_of_deg_k:Vec<usize> = g.vertices().filter(|&v| *g.phase(v).denom() == 4 && g.degree(v) == k).collect();

    'outer: for v in t_of_deg_k{
        for neigh in g.neighbors(v) {
            if *g.phase(neigh).denom() != 4 { continue 'outer;}
        }

        let mut star = vec![v]; 
        star.append(&mut g.neighbor_vec(v));
        return star;
    }

    vec![]
}


//TODO: There is a better decomposition:
fn star6_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize) -> Graph{

 //fuse out the center T , cut the wire and decompose the cat6
    match i {
        0..=2 => {
            
            let mut gn = g.clone();
            let reso0 = gn.add_vertex(VType::Z);
            let fuseout = gn.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
            gn.add_edge(reso0, fuseout);
            gn.add_to_phase(verts[0], Rational::new(-1,4));
            //normalization
            *gn.scalar_mut() *= ScalarN::sqrt2_pow(-2);
            
            //To apply the cat6
            cat_normal_form(&mut gn,&verts);
            
            cat6_replace_index(&gn, verts,i)
        },
        3..=6 => {

            let mut gn = g.clone();
            let reso1 = gn.add_vertex_with_phase(VType::Z,Rational::new(1, 1));
            let fuseout = gn.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
            gn.add_edge(reso1, fuseout);
            gn.add_to_phase(verts[0], Rational::new(1,4));

            //normalization
            *gn.scalar_mut() *= ScalarN::sqrt2_pow(-2);

            //To apply the cat6
            cat_normal_form(&mut gn,&verts);
            
            cat6_replace_index(&gn, verts,&(i-3))
        }
        other => panic!("Tried to access the index {} of the star6 decomposition",other)
    }

}


impl Decomposition {
    pub fn star6_decomp()->Decomposition{
         Decomposition { finder: |g: &Graph| -> Vec<usize>{ find_star_k(g, 6) }, nb_terms: 6, to_normal_form: trivial_normal_form, get_term: star6_replace_index, approx_alpha: |g,vert| 0.369, compute_alpha_until: THRESHOLD_COMPUTE_ALPHA }
     }
 }

 fn star3_replace_index(g: &Graph,verts: &Vec<usize>, i: &usize) -> Graph{

    //fuse out the center T , cut the wire and decompose the cat6
       match i {
           0 => {
               
               let mut gn = g.clone();
               let reso0 = gn.add_vertex(VType::Z);
               let fuseout = gn.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
               gn.add_edge(reso0, fuseout);
               gn.add_to_phase(verts[0], Rational::new(-1,4));
               //normalization
               *gn.scalar_mut() *= ScalarN::sqrt2_pow(-2);
               
               //To apply the pi copy
               cat_normal_form(&mut gn,&verts);
               gn
           },
           1 => {
   
               let mut gn = g.clone();
               let reso1 = gn.add_vertex_with_phase(VType::Z,Rational::new(1, 1));
               let fuseout = gn.add_vertex_with_phase(VType::Z,Rational::new(1, 4));
               gn.add_edge(reso1, fuseout);
               gn.add_to_phase(verts[0], Rational::new(1,4));
   
               //normalization
               *gn.scalar_mut() *= ScalarN::sqrt2_pow(-2);
   
               //To copy the pi
               cat_normal_form(&mut gn,&verts);
               
               gn
           }
           other => panic!("Tried to access the index {} of the star3 decomposition",other)
       }
   
   }

impl Decomposition {
    pub fn star3_decomp()->Decomposition{
         Decomposition { finder: |g: &Graph| -> Vec<usize>{ find_star_k(g, 3) }, nb_terms: 2, to_normal_form: trivial_normal_form, get_term: star3_replace_index, approx_alpha: |g,vert| 0.333, compute_alpha_until: THRESHOLD_COMPUTE_ALPHA }
     }
 }

