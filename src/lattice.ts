import { Lattice } from './simulation';

export default {
  dims: [],
  weights: [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36],
  opposites: [0, 2, 1, 4, 3, 6, 5, 8, 7],
  vectors: [[0,0], [0,1], [0,-1], [1,0], [-1,0], [1,1], [-1,-1], [1,-1], [-1,1]],
} as Lattice;