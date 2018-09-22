#include "advance.H"
__global__ void 
b4xsweep_kernel(pdeParam param)
{
    int col = threadIdx.x + blockIdx.x*blockDim.x;
    int row = threadIdx.y + blockIdx.y*blockDim.y;
    bool grid_valid = (row < param.cellsY && col < param.cellsX);
    if (grid_valid)
    {
        real h = param.getElement_qNew(row,col,0);	
        if (h < 0.001)
        {
            param.setElement_qNew(row, col, 0, h > 0.0 ? h : 0.0);
            param.setElement_qNew(row, col, 1, 0.0);
            param.setElement_qNew(row, col, 2, 0.0);
        }
    }
}
__global__ void 
b4ysweep_kernel(pdeParam param)
{
    int col = threadIdx.x + blockIdx.x*blockDim.x;
    int row = threadIdx.y + blockIdx.y*blockDim.y;
    bool grid_valid = (row < param.cellsY && col < param.cellsX);
    if (grid_valid)
    {
        real h = param.getElement_q_tmp(row,col,0);	
        if (h < 0.001)
        {
            param.setElement_q_tmp(row, col, 0, h > 0.0 ? h : 0.0);
            param.setElement_q_tmp(row, col, 1, 0.0);
            param.setElement_q_tmp(row, col, 2, 0.0);
        }
    }
}
