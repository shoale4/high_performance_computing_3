#include <mpi.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#define Ndim 2

/* Note: All code is provided as a potential reference serial implementation.
 *
 * You are not required to use said code
 *
 * Said code is not necessarily efficient, but is provided so that you can focus on the parallel portion
 *
 * Said code is also not necessarily representative of good C/C++ programming practices
 *
 * Said code may not necessarily lead to a nice parallel implementation
 *
 * You are able to make modifications to the code if required for your parallel implementation
 *
 * If you do find any bugs in the code, please report them on piazza
 *
 */

/** @brief Provide a default parameter for vec_to_str **/
template <typename T>
std::string vec_to_str(const std::vector<T>& vec,
                       const std::string& delim = ",");

/** @brief A function which converts a vector to a string **/
template <typename T>
std::string vec_to_str(const std::vector<T>& vec, const std::string& delim) {
  std::ostringstream oss;
  if (!vec.empty()) {
    std::copy(vec.begin(), vec.end() - 1,
              std::ostream_iterator<T>(oss, delim.c_str()));
    oss << vec.back();
  }
  return oss.str();
}

/* A sequential function to update a grid,
 * uses LDA (leading dimension) to allow for subgrid considerations */
void update_state(int m, int n, int lda, const int* in_grid, int* out_grid) {
  
  for (int i = 0; i < m; i++) { // For each row
    for (int j = 0; j < n; j++) { // For each column
      //Consider a single element
      int lin_loc = i * lda + j; //This is the linear index of the element

      int alive = 0; 

      /* Look at each neighbor */
      for (int k = -1; k < 2; k++) {
        for (int l = -1; l < 2; l++) {
          /* Figure out the index associated with each neighbor */
          int y_loc = i + k;
          int x_loc = j + l;
          int neighbor_lin_loc = y_loc * lda + x_loc;

          /* Ensure the considered neighbor is in bounds */
          if ((x_loc >= 0) && (y_loc >= 0) && (y_loc < m) && (x_loc < n)) {
            /* Check that the neighbor is actually a neighbor */
            if (!(k == 0 && l == 0)) {
              /* If it is alive, count it as alive */
              if (in_grid[neighbor_lin_loc]) {
                alive++;
              }
            }
          }
        }
      }

      /* Based on the number of alive neighbors, update the output accordingly */
      if (in_grid[lin_loc]) {
        if (alive < 2) {
          out_grid[lin_loc] = 0;
        } else if (alive > 3) {
          out_grid[lin_loc] = 0;
        } else {
          out_grid[lin_loc] = 1;
        }
      } else {
        if (alive == 3) {
          out_grid[lin_loc] = 1;
        } else {
          out_grid[lin_loc] = 0;
        }
      }
    }
  }
}

void update_state_interior(int m_loc, int n_loc, int lda_loc, int m_glb, int n_glb, const int* in_grid, int* out_grid, int coord_x, int coord_y){
  for (int i = 1; i < m_loc - 1; i++) { // For each row
    if(coord_y*m_loc + i >= m_glb){
      //std::cout << "error" << std::endl;
      continue;
    }
    for (int j = 1; j < n_loc - 1; j++) { // For each column
      //Consider a single element

      if(coord_x*n_loc + j >= n_glb){
        //std::cout << "error" << std::endl;
        continue;
      }

      int lin_loc = i * lda_loc + j; //This is the linear index of the element

      int alive = 0; 

      /* Look at each neighbor */
      for (int k = -1; k < 2; k++) {
        for (int l = -1; l < 2; l++) {
          /* Figure out the index associated with each neighbor */
          int y_loc = i + k;
          int x_loc = j + l;
          int neighbor_lin_loc = y_loc * lda_loc + x_loc;

          /* Ensure the considered neighbor is in bounds */
          if ((x_loc >= 0) && (y_loc >= 0) && (y_loc < m_loc) && (x_loc < n_loc)) {
            /* Check that the neighbor is actually a neighbor */
            if (!(k == 0 && l == 0)) {
              /* If it is alive, count it as alive */
              if (in_grid[neighbor_lin_loc]) {
                alive++;
              }
            }
          }
        }
      }

      /* Based on the number of alive neighbors, update the output accordingly */
      if (in_grid[lin_loc]) {
        if (alive < 2) {
          out_grid[lin_loc] = 0;
        } else if (alive > 3) {
          out_grid[lin_loc] = 0;
        } else {
          out_grid[lin_loc] = 1;
        }
      } else {
        if (alive == 3) {
          out_grid[lin_loc] = 1;
        } else {
          out_grid[lin_loc] = 0;
        }
      }
    }
  }
}

void update_state_boundary(int m_loc, int n_loc, int lda_loc, int m_glb, int n_glb, const int* in_grid, int* out_grid, int coord_x, int coord_y, int dim_x, int dim_y, int* halo_n, int* halo_s, int* halo_w, int* halo_e, int halo_nw, int halo_ne, int halo_sw, int halo_se){

 
  for (int i = 0; i < m_loc ; i++) {
    if(coord_y*m_loc + i >= m_glb){
      //std::cout << "error" << std::endl;
      continue;
    }
    for (int j = 0; j < n_loc ; j++) { // For each column
      //Consider a single element

      if(coord_x*n_loc + j >= n_glb){
        //std::cout << "error" << std::endl;
        continue;
      }

      if(i > 0 && i < m_loc - 1 && j > 0 && j < n_loc - 1){
        continue;
      }

      int lin_loc = i * lda_loc + j; //This is the linear index of the element

      int alive = 0; 

      /* Look at each neighbor */
      for (int k = -1; k < 2; k++) {
        for (int l = -1; l < 2; l++) {
          /* Figure out the index associated with each neighbor */
          int y_loc = i + k;
          int x_loc = j + l;
          int neighbor_lin_loc = y_loc * lda_loc + x_loc;


          if(y_loc == -1 && x_loc == -1){
            if(coord_x > 0 && coord_y > 0){
              if(halo_nw){
                alive++;
              }              
            }
          }
          else if(y_loc == -1 && x_loc == n_loc ){
            if(coord_x < dim_x - 1 && coord_y > 0){
              if(halo_ne){
                alive++;
              }
            }
          }
          else if(y_loc == m_loc && x_loc == -1 ){
            if(coord_x > 0 && coord_y < dim_y - 1){
              if(halo_sw){
                alive++;
              }
            }
          }
          else if(y_loc == m_loc && x_loc == n_loc ){
            if(coord_x < dim_x - 1 && coord_y < dim_y - 1){
              if(halo_se){
                alive++;
              }             
            }

          }
          else if(y_loc == -1){
            if(coord_y > 0){
              if(halo_n[x_loc]){
                alive++;
              }
            }
          }
          else if(y_loc == m_loc){
            if(coord_y < dim_y - 1){
              if(halo_s[x_loc]){
                alive++;
              }
            }
          }
          else if(x_loc == -1){
            if(coord_x > 0){
              if(halo_w[y_loc]){
                alive++;
              }
            }
          }
          else if(x_loc == n_loc){
            if(coord_x < dim_x - 1){
              if(halo_e[y_loc]){
                alive++;
              }
            }
          }
          else{
              if (!(k == 0 && l == 0)) {
                if (in_grid[neighbor_lin_loc]) {
                  alive++;
                }
            }
          }
        }
      }

      /* Based on the number of alive neighbors, update the output accordingly */
      if (in_grid[lin_loc]) {
        if (alive < 2) {
          out_grid[lin_loc] = 0;
        } else if (alive > 3) {
          out_grid[lin_loc] = 0;
        } else {
          out_grid[lin_loc] = 1;
        }
      } else {
        if (alive == 3) {
          out_grid[lin_loc] = 1;
        } else {
          out_grid[lin_loc] = 0;
        }
      }
    }
  }
}

/* Read in the data from an input file */
void read_data(const std::string& input_filename, int m, int n,
               std::vector<int>& output_data) {
  output_data.reserve(m * n);
  std::ifstream input_file(input_filename, std::ios::in);
  for (int i = 0; i < m * n; i++) {
    std::string mystring;
    input_file >> mystring;
    output_data.push_back(std::stoi(mystring));
  }
}

/* Write output data to a file */
void write_data(const std::string& output_filename, int m, int n,
                const std::vector<int>& output_data) {
  std::ofstream output_file(output_filename, std::ios::out);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      output_file << output_data[i * n + j];
      if (j < n - 1) {
        output_file << " ";
      }
    }
    output_file << "\n";
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  const int m = std::stoi(argv[1]);
  const int n = std::stoi(argv[2]);
  const int gen = std::stoi(argv[3]);
  const std::string input_file(argv[4]);
  const std::string output_file(argv[5]);

  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> global_data;

  /* Allocate the output data buffer */
  std::vector<int> output_data;








  if (size == 1) { /* If serial, use the serial code */

    /* On root, read in the data */
    if (rank == 0) {
      read_data(input_file, m, n, global_data);
    }



    output_data.reserve(global_data.size());

    /* For each generation update the state */
    for (int i = 0; i < gen; i++) {
      update_state(m, n, n, global_data.data(), output_data.data());

      /* Swap the input and output */
      if (i < gen - 1) {
        std::swap(global_data, output_data);
      }
    }

    /* On root, output the data */
    if (rank == 0) {
      write_data(output_file, m, n, output_data);
    }
  } 

  else {
    /* You implement this */

    int dims[Ndim] = {0, 0};
    int period[Ndim] = {0, 0};
    MPI_Comm grid_comm;
    int grid_rank;
    int grid_coords[Ndim];
    int xdim, ydim;
    int m_pad, n_pad;


    MPI_Dims_create(size, Ndim, dims);
    MPI_Cart_create(MPI_COMM_WORLD, Ndim, dims, period, 1, &grid_comm);
    MPI_Comm_rank(grid_comm, &grid_rank);
    MPI_Cart_coords(grid_comm, grid_rank, Ndim, grid_coords);

    //std::cout << grid_rank << grid_coords[0] << grid_coords[1] << std::endl;

    if(grid_rank == 0){
      read_data(input_file, m, n, global_data);
    }


    n_pad = ((ceil((float)(n*1.0/dims[0]))))*dims[0];
    m_pad = ((ceil((float)(m*1.0/dims[1]))))*dims[1];

    //std::cout << m_pad << n_pad << std::endl;

    xdim = (int)(n_pad/dims[0]);
    ydim = (int)(m_pad/dims[1]);

    std::vector<int> global_data_pad(m_pad*n_pad,0);
    std::vector<int> output_data_pad;

    if(grid_rank == 0){
      for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
          global_data_pad[i*n_pad + j] = global_data[i*n + j];
        }
      }
    }

    output_data.reserve(global_data.size());
    output_data_pad.reserve(global_data_pad.size());



   // std::cout << xdim << " " << ydim << std::endl;

    std::vector<int> local_data(ydim*xdim);
    std::vector<int> local_data_out; 
    local_data_out.reserve(local_data.size());

    int counts[dims[1]*dims[0]];
    int displs[dims[1]*dims[0]];

    for(int j = 0; j < dims[0]; j++){
      for(int i = 0; i < dims[1]; i++){
        counts[j*dims[1] + i] = 1;
        displs[j*dims[1] + i] = i*xdim*ydim*dims[0] + j*xdim;
        //displs[j*dims[1] + i] = 8;
      }
    }

    //std::cout << displs[0] << displs[1] << displs[2] << displs[3] << std::endl;



    MPI_Datatype block, temp_block;
    MPI_Type_vector(ydim, xdim, n_pad, MPI_INT, &temp_block);
    MPI_Type_create_resized(temp_block, 0, sizeof(int), &block);
    MPI_Type_commit(&block);


    MPI_Scatterv(&global_data_pad[0], counts, displs, block, &local_data[0], ydim*xdim, MPI_INT, 0, grid_comm);



    MPI_Comm row_comm, col_comm;
    int remain_dims_row[Ndim] = {false, true};
    int remain_dims_col[Ndim] = {true, false};
    MPI_Cart_sub(grid_comm, remain_dims_row, &row_comm);
    MPI_Cart_sub(grid_comm, remain_dims_col, &col_comm);
    int row_rank, row_size, col_rank, col_size;
    MPI_Comm_size(row_comm, &row_size);
    MPI_Comm_size(col_comm, &col_size);
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);



    int halo_n[xdim], halo_s[xdim], halo_w[ydim], halo_e[ydim];
    int halo_nw, halo_ne, halo_sw, halo_se ;

    MPI_Datatype col, temp_col;
    MPI_Type_vector(ydim, 1, xdim, MPI_INT, &temp_col);
    MPI_Type_create_resized(temp_col, 0, sizeof(int), &col);
    MPI_Type_commit(&col);

    int N, W, S, E, NW = -1, SE = -1, SW = -1, NE = -1;
    N = row_rank - 1;
    S = row_rank + 1;
    W = col_rank - 1;
    E = col_rank + 1;

    int coords_NW[2] = {grid_coords[0] - 1, grid_coords[1] - 1};
    int coords_NE[2] = {grid_coords[0] + 1, grid_coords[1] - 1};
    int coords_SW[2] = {grid_coords[0] - 1, grid_coords[1] + 1};
    int coords_SE[2] = {grid_coords[0] + 1, grid_coords[1] + 1};

    if(grid_coords[0] > 0 && grid_coords[1] > 0){
      MPI_Cart_rank(grid_comm, coords_NW, &NW);     
    }

    if(grid_coords[0] < dims[0] - 1 && grid_coords[1] > 0){
      MPI_Cart_rank(grid_comm, coords_NE, &NE);
    }

    if(grid_coords[0] > 0 && grid_coords[1] < dims[1] - 1){
      MPI_Cart_rank(grid_comm, coords_SW, &SW);
    }
    if(grid_coords[0] < dims[0] - 1 && grid_coords[1] < dims[1] - 1){
      MPI_Cart_rank(grid_comm, coords_SE, &SE); 
    }

    MPI_Request send_req_n, send_req_s, send_req_w, send_req_e, send_req_nw, send_req_ne, send_req_sw, send_req_se;
    MPI_Request recv_req_n, recv_req_s, recv_req_w, recv_req_e, recv_req_nw, recv_req_ne, recv_req_sw, recv_req_se;

    
    if(row_rank > 0){
      MPI_Send_init(&local_data[0], xdim, MPI_INT, N, 11, row_comm, &send_req_n);
      MPI_Recv_init(&halo_n,xdim, MPI_INT, N, 11, row_comm, &recv_req_n);
    }

    if(row_rank < dims[1] - 1){
      MPI_Send_init(&local_data[(ydim - 1)*xdim], xdim, MPI_INT, S, 11, row_comm, &send_req_s);
      MPI_Recv_init(&halo_s,xdim, MPI_INT, S, 11, row_comm, &recv_req_s);
    }
    
    if(col_rank > 0){
      MPI_Send_init(&local_data[0], 1, col, W, 11, col_comm, &send_req_w);
      MPI_Recv_init(&halo_w,ydim, MPI_INT, W, 11, col_comm, &recv_req_w);
    }
    
    if(col_rank < dims[0] - 1){
      MPI_Send_init(&local_data[xdim - 1], 1, col, E, 11, col_comm, &send_req_e);
      MPI_Recv_init(&halo_e,ydim, MPI_INT, E, 11, col_comm, &recv_req_e);
    }


    if(grid_coords[0] > 0 && grid_coords[1] > 0){
       MPI_Send_init(&local_data[0], 1, MPI_INT, NW, 11, grid_comm, &send_req_nw); 
       MPI_Recv_init(&halo_nw, 1, MPI_INT, NW, 11, grid_comm, &recv_req_nw);   
    }

    if(grid_coords[0] < dims[0] - 1 && grid_coords[1] > 0){
      MPI_Send_init(&local_data[xdim - 1], 1, MPI_INT, NE, 11, grid_comm, &send_req_ne);
      MPI_Recv_init(&halo_ne, 1, MPI_INT, NE, 11, grid_comm, &recv_req_ne);
    }

    if(grid_coords[0] > 0 && grid_coords[1] < dims[1] - 1){
      MPI_Send_init(&local_data[(ydim - 1)*xdim], 1, MPI_INT, SW, 11, grid_comm, &send_req_sw);
      MPI_Recv_init(&halo_sw, 1, MPI_INT, SW, 11, grid_comm, &recv_req_sw);

    }
    if(grid_coords[0] < dims[0] - 1 && grid_coords[1] < dims[1] - 1){
      MPI_Send_init(&local_data[ydim*xdim - 1], 1, MPI_INT, SE, 11, grid_comm, &send_req_se);
      MPI_Recv_init(&halo_se, 1, MPI_INT, SE, 11, grid_comm, &recv_req_se);
    }

    for (int i = 0; i < gen; i++) {

      if(row_rank > 0){
        MPI_Start(&send_req_n);
        MPI_Start(&recv_req_n);
      }

      if(row_rank < dims[1] - 1){
        MPI_Start(&send_req_s);
        MPI_Start(&recv_req_s);    }
      
      if(col_rank > 0){
        MPI_Start(&send_req_w);
        MPI_Start(&recv_req_w);
      }
      
      if(col_rank < dims[0] - 1){
        MPI_Start(&send_req_e);
        MPI_Start(&recv_req_e);
      }


      if(grid_coords[0] > 0 && grid_coords[1] > 0){
        MPI_Start(&send_req_nw);
        MPI_Start(&recv_req_nw); 
      }

      if(grid_coords[0] < dims[0] - 1 && grid_coords[1] > 0){
        MPI_Start(&send_req_ne);
        MPI_Start(&recv_req_ne);
      }

      if(grid_coords[0] > 0 && grid_coords[1] < dims[1] - 1){
        MPI_Start(&send_req_sw);
        MPI_Start(&recv_req_sw);

      }
      if(grid_coords[0] < dims[0] - 1 && grid_coords[1] < dims[1] - 1){
        MPI_Start(&send_req_se);
        MPI_Start(&recv_req_se);
      }

      update_state_interior(ydim, xdim, xdim, m, n, local_data.data(), local_data_out.data(), grid_coords[0], grid_coords[1]);

      if(row_rank > 0){
        MPI_Wait(&send_req_n, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_n, MPI_STATUS_IGNORE);
      }

      if(row_rank < dims[1] - 1){
        MPI_Wait(&send_req_s, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_s, MPI_STATUS_IGNORE);    }
      
      if(col_rank > 0){
        MPI_Wait(&send_req_w, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_w, MPI_STATUS_IGNORE);
      }
      
      if(col_rank < dims[0] - 1){
        MPI_Wait(&send_req_e, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_e, MPI_STATUS_IGNORE);
      }


      if(grid_coords[0] > 0 && grid_coords[1] > 0){
        MPI_Wait(&send_req_nw, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_nw, MPI_STATUS_IGNORE); 
      }

      if(grid_coords[0] < dims[0] - 1 && grid_coords[1] > 0){
        MPI_Wait(&send_req_ne, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_ne, MPI_STATUS_IGNORE);
      }

      if(grid_coords[0] > 0 && grid_coords[1] < dims[1] - 1){
        MPI_Wait(&send_req_sw, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_sw, MPI_STATUS_IGNORE);

      }
      if(grid_coords[0] < dims[0] - 1 && grid_coords[1] < dims[1] - 1){
        MPI_Wait(&send_req_se, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_req_se, MPI_STATUS_IGNORE);
      }

/*      if(grid_rank == 1){
        std::cout << halo_ne << std::endl;
        std::cout << halo_n[0] << " " << halo_n[1] << std::endl;
        std::cout << halo_e[0] << " " << halo_e[1] << std::endl;
      }
*/
      update_state_boundary(ydim, xdim, xdim, m, n, local_data.data(), local_data_out.data(), grid_coords[0], grid_coords[1], dims[0], dims[1], halo_n, halo_s, halo_w, halo_e, halo_nw, halo_ne, halo_sw, halo_se);

      
      if (i < gen - 1) {
        //std::swap(local_data, local_data_out);
        for(int j = 0; j < ydim; j++){
          for(int k = 0; k < xdim; k++){
            local_data[j*xdim + k] = local_data_out[j*xdim + k];
          }
        }
      }

      MPI_Barrier(grid_comm);



    }

    
  

    MPI_Gatherv(&local_data_out[0], ydim*xdim, MPI_INT, &output_data_pad[0], &counts[0], &displs[0], block, 0, grid_comm);    



    if(grid_rank == 0){
      for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
          output_data[i*n + j] = output_data_pad[i*n_pad + j];
        }
      }
    }


    if (grid_rank == 0) {
      write_data(output_file, m, n, output_data);
    }  
 
  }

  MPI_Finalize();
}
