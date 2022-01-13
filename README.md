# FiniteDifference
This code uses finite difference and forward euler methods to solve allen-cahn differential equations.


main.py: Use this code to compute the trajectory of the vector field. The code is used to set initial parameters and run the methods in methods.py.
  

Methods.py: Contains the FiniteDifference and SquaredMesh classes.

  FiniteDifference(method = 'central')
     Uses finite difference method to evaluate the evolution of a function subjected to the Allen-Cahn equation dC / dt = -L grad(F'), where C is the concentration vector field, F is the free energy functional, and L is a constant.
    
   Attribuites:
            method: 'forward', 'backward', or 'central'   

   Methods:
            
            
            backward(vector,dn, order = 1):
                uses the forward finite difference to compute first and second order derivatives using df(x)/dx = f(x) - f(x-h) / delta_x and for first order derivatives and df(x)/dx = f(x-2h) + 2*f(x) - f(x - h) / delta_x ˆ2 for second order derivatives, and using dirichlet boundary conditions.
                - vector: array float type. Must contain the values for every point in the mesh grid for the problem.
                - dn: float type. The finite position difference between points in the mesh grid.
                - order: int type, Default = 1. The derivative order outputed from the function. Must be 1 or two for the current version
                
            forward(vector,dn, order = 1): 
                uses the forward finite difference to compute first and second order derivatives using dF(x)/dx = f(x+h) - f(x) / delta_x and for first order derivatives and df(x)/dx = f(x+2h) + 2*f(x) - f(x + h) / delta_x ˆ2 for second order derivatives, and using dirichlet boundary conditions.
                - vector: array float type. Must contain the values for every point in the mesh grid for the problem.
                - dn: float type. The finite position difference between points in the mesh grid.
                - order: int type, Default = 1. The derivative order outputed from the function. Must be 1 or two for the current version
            
            central(vector,dn, order = 1): 
                uses the forward finite difference to compute first and second order derivatives using df(x)/dx = f(x+h) - f(x -h) / delta_x for first order derivatives and df(x)/dx = f(x+h) - 2*f(x) + f(x - h) / delta_x ˆ2 for second order derivatives and using dirichlet boundary conditions.
                - vector: array float type. Must contain the values for every point in the mesh grid for the problem.
                - dn: float type. The finite position difference between points in the mesh grid.
                - order: int type, Default = 1. The derivative order outputed from the function. Must be 1 or two for the current version
            
            onedimensional(c, dt, dx, L = 1, alpha = 1, beta = 1, Nt = 50): 
                Uses forward euler method to integrate over time a one dimensional vector C
                - c: array float type. Contains the values for c in every mesh point.
                - dt: float type. This is the time step size.
                - dx: float type. TThe finite position difference between points in the mesh grid in the x direction. 
                - L: float type, Default = 1. This is a constant in the allen-cahn differential equation.
                - alpha: float type, Default = 1. This is a constant in the free energy functional.
                - beta: float type, Default = 1. This is a constant in the free energy functional. 
                - Nt: int type, Default = 50. The number of time steps to be performed.
            
            twodimensional(c, dt, dx, dy, L = 1, alpha = 1, beta = 1, Nt = 50): 
                Uses forward euler method to integrate over time a two dimensional vector C
                - c: array float type. Contains the values for c in every mesh point.
                - dt: float type. This is the time step size.
                - dx: float type. TThe finite position difference between points in the mesh grid in the x direction. 
                - dy: float type. TThe finite position difference between points in the mesh grid in the y direction. 
                - L: float type, Default = 1. This is a constant in the allen-cahn differential equation.
                - alpha: float type, Default = 1. This is a constant in the free energy functional.
                - beta: float type, Default = 1. This is a constant in the free energy functional. 
                - Nt: int type, Default = 50. The number of time steps to be performed.
            
            threedimensional(c, dt, dx, dy, dz, L = 1, alpha = 1, beta = 1, Nt = 50): 
                Uses forward euler method to integrate over time a three dimensional vector C
                - c: array float type. Contains the values for c in every mesh point.
                - dt: float type. This is the time step size.
                - dx: float type. TThe finite position difference between points in the mesh grid in the x direction. 
                - dy: float type. TThe finite position difference between points in the mesh grid in the y direction. 
                - dz: float type. TThe finite position difference between points in the mesh grid in the z direction. 
                - L: float type, Default = 1. This is a constant in the allen-cahn differential equation.
                - alpha: float type, Default = 1. This is a constant in the free energy functional.
                - beta: float type, Default = 1. This is a constant in the free energy functional. 
                - Nt: int type, Default = 50. The number of time steps to be performed.
