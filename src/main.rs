use noisy_float::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs;
extern crate tsplib;
use tsplib::*;

fn main() {
  println!("Hello, TSP. Using testdata/berlin52.tsp as input !");
  // let mut input = String::new();
  // io::stdin().read_to_string(&mut input);
  let mut tsp = TSP::new("testdata/berlin52.tsp".to_string());
  tsp.solve();
}
// simple exact TSP solver based on branch-and-bound/Held--Karp
#[derive(Debug, Clone)]
struct TSP {
  n: usize,
  x: Vec<N32>,
  y: Vec<N32>,
  cost: Vec<Vec<N32>>,
  cost_with_pi: Vec<Vec<N32>>,
  best: Node,
}
#[derive(Eq, Default, Debug, Clone)]
struct Node {
  excluded: Vec<Vec<bool>>,
  pi: Vec<N32>,
  lower_bound: N32,
  degree: Vec<usize>,
  parent: Vec<usize>,
}
impl Ord for Node {
  fn cmp(&self, other: &Node) -> Ordering {
    self.lower_bound.cmp(&other.lower_bound)
  }
}
impl PartialEq for Node {
  fn eq(&self, other: &Self) -> bool {
    self.lower_bound == other.lower_bound
  }
}
impl PartialOrd for Node {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}
impl TSP {
  fn add_edge(&self, node: &mut Node, i: usize, j: usize) {
    node.lower_bound += self.cost_with_pi[i][j];
    node.degree[i] += 1;
    node.degree[j] += 1;
  }
  fn exclude(&mut self, node: &mut Node, i: usize, j: usize) -> Node {
    let mut child: Node = Default::default();
    child.pi = vec![n32(0.0); self.n];
    child.parent = vec![0; self.n];
    child.excluded = node.excluded.clone();
    child.excluded[i] = node.excluded[i].clone();
    child.excluded[j] = node.excluded[j].clone();
    child.excluded[i][j] = true;
    child.excluded[j][i] = true;
    self.compute_held_karp(&mut child);
    child
  }
  fn compute_held_karp(&mut self, node: &mut Node) {
    node.lower_bound = n32(std::f32::MIN);
    node.degree = vec![0; self.n];
    node.parent = vec![0; self.n];
    let mut lambda = n32(0.1);
    while lambda > 1e-06 {
      let previous_lower = node.lower_bound;
      self.compute_one_tree(node);
      if node.lower_bound >= self.best.lower_bound {
        return;
      }
      if node.lower_bound >= previous_lower {
        lambda *= 0.9;
      }
      let mut denom = 0;
      for i in 1..self.n {
        let deg2 = (node.degree[i as usize] as i32) - 2;
        denom += deg2 * deg2;
      }
      // println!("---DENOM : {}", denom);
      if denom == 0 {
        return;
      }
      // println!("---lambda, lower_bound : {}, {}", lambda, node.lower_bound);
      let t: N32 = lambda * node.lower_bound / n32(denom as f32);
      for i in 1..self.n {
        node.pi[i as usize] += t * n32(((node.degree[i as usize] as i32) - 2) as f32);
      }
    }
  }
  fn compute_one_tree(&mut self, node: &mut Node) {
    node.lower_bound = n32(0.0);
    node.degree = vec![0; self.n];
    for i in 0..self.n {
      for j in 0..self.n {
        self.cost_with_pi[i][j] = if node.excluded[i][j] {
          n32(std::f32::MAX)
        } else {
          self.cost[i][j] + node.pi[i] + node.pi[j]
        }
      }
    }
    // find the two cheapest edges from 0
    let (mut first_neighbor, mut second_neighbor) =
      if self.cost_with_pi[0][2] < self.cost_with_pi[0][1] {
        (2, 1)
      } else {
        (1, 2)
      };
    //find the top two smallest edges from 0, keeping track of the cheapest and 2nd cheapest.
    for j in 3..self.n {
      let j_cost = self.cost_with_pi[0][j];
      if j_cost < self.cost_with_pi[0][second_neighbor] {
        if j_cost < self.cost_with_pi[0][first_neighbor] {
          second_neighbor = first_neighbor;
          first_neighbor = j;
        } else {
          second_neighbor = j;
        }
      }
    }
    self.add_edge(node, 0, first_neighbor);
    node.parent = vec![first_neighbor; self.n];
    node.parent[first_neighbor] = 0;
    // compute the minimum spanning tree on nodes 1..n-1
    let mut min_cost = self.cost_with_pi[first_neighbor].clone();
    for _k in 2..self.n {
      let mut i = node.degree.iter().position(|&degree| degree == 0).unwrap();
      for j in (i + 1)..self.n {
        if node.degree[j] == 0 && min_cost[j] < min_cost[i] {
          i = j;
        }
      }
      self.add_edge(node, node.parent[i], i);
      for j in 1..self.n {
        if node.degree[j] == 0 && self.cost_with_pi[i][j] < min_cost[j] {
          min_cost[j] = self.cost_with_pi[i][j];
          node.parent[j] = i;
        }
      }
    }
    self.add_edge(node, 0, second_neighbor);
    node.parent[0] = second_neighbor;
    // println!(
    //   "Setting lower bound from {} to {}",
    //   node.lower_bound,
    //   node.lower_bound.round(),
    // );
    node.lower_bound = node.lower_bound.round();
  }
  fn solve(&mut self) {
    self.best = self.new_node();
    self.best.lower_bound = n32(std::f32::MAX);
    let mut current_node: Node = self.new_node();
    println!("{:?}", self);
    self.compute_held_karp(&mut current_node);
    let mut pq = BinaryHeap::new();
    loop {
      loop {
        let mut iopt: Option<usize> = None;
        for j in 0..self.n {
          if current_node.degree[j] > 2
            && (iopt.is_none() || current_node.degree[j] < current_node.degree[iopt.unwrap()])
          {
            iopt = Some(j);
          }
        }
        match iopt {
          None => {
            if current_node.lower_bound < self.best.lower_bound {
              self.best = current_node.clone();
              // println!("{}", self.best_lb)
            }
            break;
          }
          Some(i) => {
            println!(".");
            let mut children: BinaryHeap<Node> = BinaryHeap::new();
            let parent_i = current_node.parent[i];
            children.push(self.exclude(&mut current_node, i, parent_i));
            for j in 0..self.n {
              if current_node.parent[j] == i {
                children.push(self.exclude(&mut current_node, i, j));
              }
            }
            current_node = children.pop().unwrap();
            pq.append(&mut children);
            if current_node.lower_bound >= self.best.lower_bound {
              break;
            }
          }
        }
      }
      match pq.pop() {
        None => {
          println!("Breaking because pq.pop returned None");
          break;
        }
        Some(new_node) => current_node = new_node,
      };
      // println!("There are {} nodes to visit", pq.len());
      if current_node.lower_bound > self.best.lower_bound {
        println!(
          "Breaking because current node lower bound {} is greater than best lower bound of {}",
          current_node.lower_bound, self.best.lower_bound
        );
        break;
      }
    }
    println!("best lower bound is {}", self.best.lower_bound);
    let mut j = 0;
    let mut i = self.best.parent[0];
    let mut cost = n32(0.0);
    while i != 0 {
      cost += self.cost_with_pi[i][j];
      i = self.best.parent[j];
      println!(
        "{}->{}\t{}\t{}\t{}\t{}",
        j,
        i,
        self.x[j],
        self.y[j],
        self.x[i] - self.x[j],
        self.y[i] - self.y[j]
      );
      j = i;
    }
    println!("The trip cost is {}", cost)
  }
  fn new(input: String) -> TSP {
    let file_contents = fs::read_to_string(input).unwrap();
    let problem = parse_whole_problem(&file_contents).unwrap().1;
    let n: usize = problem.header.dimension as usize;

    let mut cost = vec![vec![n32(0.0); n]; n];
    let x: Vec<N32> = problem
      .clone()
      .data
      .node_coordinates
      .unwrap()
      .iter()
      .map(|c| c.1.clone())
      .collect();
    let y: Vec<N32> = problem
      .clone()
      .data
      .node_coordinates
      .unwrap()
      .iter()
      .map(|c| c.2.clone())
      .collect();
    // TSPLIB distances are rounded to the nearest integer to avoid the sum of square roots problem
    for i in 0..n {
      for j in 0..n {
        let dx = x[i] - x[j];
        let dy = y[i] - y[j];
        cost[i][j] = (dx * dx + dy * dy).sqrt().trunc();
      }
    }
    TSP {
      n,
      x,
      y,
      cost,
      cost_with_pi: vec![vec![n32(0.0); n]; n],
      best: Default::default(),
    }
  }
  fn new_node(&self) -> Node {
    Node {
      excluded: vec![vec![false; self.n]; self.n],
      pi: vec![n32(0.0); self.n],
      lower_bound: n32(std::f32::MAX),
      degree: vec![0; self.n],
      parent: vec![0; self.n],
    }
  }
}
