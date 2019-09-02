use noisy_float::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::env;
use std::fs;
extern crate tsplib;
use tsplib::*;

fn main() {
  let args: Vec<String> = env::args().collect();
  let name: &str = if args.len() < 2 {
    &"tests/testdata/berlin52.tsp"
  } else {
    &args[1]
  };
  println!("Using '{:}' as input !", name);
  let mut tsp = TSP::new(name);
  tsp.solve();
}
// simple exact TSP solver based on branch-and-bound/Held--Karp
#[derive(Debug, Clone)]
struct TSP {
  n: usize,
  x: Vec<N64>,
  y: Vec<N64>,
  cost: Vec<Vec<N64>>,
  cost_with_pi: Vec<Vec<N64>>,
  best: Node,
}
#[derive(Eq, Default, Debug, Clone)]
struct Node {
  excluded: Vec<Vec<bool>>,
  pi: Vec<N64>,
  lower_bound: N64,
  degree: Vec<usize>,
  parent: Vec<usize>,
}
impl Node {
  fn new(n: usize) -> Node {
    Node {
      excluded: vec![vec![false; n]; n],
      pi: vec![n64(0.0); n],
      lower_bound: n64(std::f64::MAX),
      degree: vec![0; n],
      parent: vec![0; n],
    }
  }
}
impl Ord for Node {
  fn cmp(&self, other: &Node) -> Ordering {
    self.lower_bound.cmp(&other.lower_bound).reverse()
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
    let mut child: Node = Node::new(self.n);
    child.pi = vec![n64(0.0); self.n];
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
    node.lower_bound = n64(std::f64::MIN_POSITIVE);
    node.degree = vec![0; self.n];
    node.parent = vec![0; self.n];
    let mut lambda = n64(0.1);
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
        let deg2 = (node.degree[i] as i64) - 2;
        denom += deg2 * deg2;
      }
      // println!("---DENOM : {}", denom);
      if denom == 0 {
        return;
      }
      // println!("---lambda, lower_bound : {}, {}", lambda, node.lower_bound);
      let t: N64 = lambda * node.lower_bound / n64(denom as f64);
      for i in 1..self.n {
        node.pi[i] += t * n64(((node.degree[i] as i64) - 2) as f64);
      }
    }
  }
  fn compute_one_tree(&mut self, node: &mut Node) {
    node.lower_bound = n64(0.0);
    node.degree = vec![0; self.n];
    for i in 0..self.n {
      for j in 0..self.n {
        self.cost_with_pi[i][j] = if node.excluded[i][j] {
          n64(std::f64::MAX)
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
      //Find the first degree = 0 node
      let mut i = 1;
      while i < self.n {
        if node.degree[i] == 0 {
          break;
        }
        i += 1;
      }
      // node.degree .iter() .skip(1) .position(|&degree| degree == 0) .unwrap();
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
    node.lower_bound = (node.lower_bound + 0.5).trunc();
  }
  fn solve(&mut self) {
    self.best = Node::new(self.n);
    // self.best.lower_bound = n64(std::f64::MAX);
    let mut current_node: Node = Node::new(self.n);
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
        if iopt.is_none() {
          if current_node.lower_bound < self.best.lower_bound {
            self.best = current_node.clone();
            print!("{:}", self.best.lower_bound);
            // println!(
            //   "Updating self.best because current.lower {:} < best.lower {:}",
            //   current_node.lower_bound, self.best.lower_bound
            // );
          }
          break;
        }
        print!(".");
        let i = iopt.unwrap();
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
        // println!(
        //       "Queue length is {:}, current lower bound is {:}, best lower bound is {:}, total degree is {:}",
        //       pq.len(),
        //       current_node.lower_bound,
        //       self.best.lower_bound,
        //       current_node.degree.iter().sum::<usize>()
        //     );
        if current_node.lower_bound >= self.best.lower_bound {
          break;
        }
      }
      println!("");
      match pq.pop() {
        None => {
          println!("Breaking because the solution heap is empty");
          break;
        }
        Some(new_node) => current_node = new_node,
      };
      if current_node.lower_bound > self.best.lower_bound {
        println!(
          "Breaking because current node lower bound {} is greater than best lower bound of {}",
          current_node.lower_bound, self.best.lower_bound
        );
        break;
      }
    }
    self.print_sol(&self.best);
  }

  fn print_sol(&self, node: &Node) {
    println!("best lower bound is {}", node.lower_bound);
    let mut j = 0;
    let mut i = node.parent[0];
    let mut cost = n64(0.0);
    while i != 0 {
      i = node.parent[j];
      println!(
        "{}->{}\t{}\t{}\t{}\t{}\t{}",
        j + 1,
        i + 1,
        self.x[j],
        self.y[j],
        self.x[i] - self.x[j],
        self.y[i] - self.y[j],
        self.cost[i][j],
      );
      cost += self.cost[i][j];
      j = i;
    }
    println!("The trip cost is {}", cost)
  }

  fn new(input: &str) -> TSP {
    let file_contents = fs::read_to_string(input).unwrap();
    let problem = parse_whole_problem(&file_contents).unwrap().1;
    let n: usize = problem.header.dimension as usize;

    let mut cost = vec![vec![n64(0.0); n]; n];
    let x: Vec<N64> = problem
      .clone()
      .data
      .node_coordinates
      .unwrap()
      .iter()
      .map(|c| n64(c.1.raw().into()))
      .collect();
    let y: Vec<N64> = problem
      .clone()
      .data
      .node_coordinates
      .unwrap()
      .iter()
      .map(|c| n64(c.2.raw().into()))
      .collect();
    // TSPLIB distances are rounded to the nearest integer to avoid the sum of square roots problem
    for i in 0..n {
      for j in 0..n {
        let dx = x[i] - x[j];
        let dy = y[i] - y[j];
        cost[i][j] = ((dx * dx + dy * dy).sqrt() + 0.5).trunc();
      }
    }
    TSP {
      n,
      x,
      y,
      cost,
      cost_with_pi: vec![vec![n64(0.0); n]; n],
      best: Node::new(n),
    }
  }
}
