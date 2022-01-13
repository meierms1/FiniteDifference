#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 14:49:05 2022

@author: maycon
"""
import numpy as np


''' 
FiniteDifference: Uses finite difference method to evaluate the evolution of a function 
            subjected to the Allen-Cahn equation dC / dt = -L grad(F'), where C is the 
            concentration vector field, F is the free energy functional, and L is a constant.

        Attribuites
            method: 'forward', 'backward', or 'central'   

        Methods
            forward: uses the forward finite difference to compute first and 
                    second order derivatives using dF(x)/dx = f(x+h) - f(x) / delta_x and 
                    for first order derivatives and df(x)/dx = f(x+2h) + 2*f(x) - f(x + h) / delta_x ˆ2
                    for second order derivatives, and using dirichlet boundary conditions
            backward:uses the forward finite difference to compute first and 
                    second order derivatives using df(x)/dx = f(x) - f(x-h) / delta_x and for
                    first order derivatives and df(x)/dx = f(x-2h) + 2*f(x) - f(x - h) / delta_x ˆ2
                    for second order derivatives, and using dirichlet boundary conditions
            central: uses the forward finite difference to compute first and 
                    second order derivatives using df(x)/dx = f(x+h) - f(x -h) / delta_x for
                    first order derivatives and df(x)/dx = f(x+h) - 2*f(x) + f(x - h) / delta_x ˆ2
                    for second order derivatives and using dirichlet boundary conditions.
            onedimensional: Uses forward euler method to integrate over time a one dimensional vector C
            twodimensional: Uses forward euler method to integrate over time a two dimensional vector C
            threedimensional: Uses forward euler method to integrate over time a three dimensional vector C                                                 

'''

class FiniteDifference:
    def __init__(self, method = 'central'):        
        self.method = method
        
    '''
    The forward, backward and central fucntions use the dirichlet boundary condition
    '''
    def forward(self, vector,dn, order = 1):
        self.order = order
        number_of_elements = len(vector)
        dfdn = np.zeros(number_of_elements)
        if self.order == 1:
            for i in range(number_of_elements):
                if i == number_of_elements:
                    f = vector[i]
                    f_forward = 0
                else:
                    f = vector[i]
                    f_forward = vector[i+1] 
                dfdn[i] = (f_forward - f)/ dn
            return dfdn
        
        elif self.order == 2:
            for i in range(number_of_elements):
                if i == number_of_elements:
                    f = vector[i]
                    f_forward = 0
                    f_forward2 = 0
                elif i == number_of_elements - 1:
                    f = vector[i]
                    f_forward = vector[i+1]
                    f_forward2 = 0    
                else:
                    f = vector[i]
                    f_forward = vector[i+1]
                    f_forward2 = vector[i+2]
                dfdn[i] = (f_forward2 - 2 * f_forward + f)/ (dn ** 2)                
            return dfdn
        else:
            raise ValueError('this order derivative not implemented')
        
    def backward(self, vector, dn, order = 1):
        self.order = order
        number_of_elements = len(vector)
        dfdn = np.zeros(number_of_elements)
        if self.order == 1:
            for i in range(number_of_elements):
                if i == 0:
                    f = vector[i]
                    f_back = 0
                else:
                    f = vector[i]
                    f_back = vector[i-1]
                dfdn[i] = (f - f_back) / dn
            return dfdn
        
        elif self.order == 2:
            for i in range(number_of_elements):
                if i == 0:
                    f = vector[i]
                    f_back = 0
                    f_back2 = 0
                elif i == 1:
                    f = vector[i]
                    f_back = vector[i-1]
                    f_back2 = 0
                else:
                    f = vector[i]
                    f_back = vector[i-1]
                    f_back2 = vector[i-2]
                    
                dfdn[i] = (f + f_back2 - 2 * f_back) / (dn ** 2)
                
            return dfdn
        else:
            raise ValueError('this order derivative is not implemented')
    
    def central(self, vector, dn, order = 1):
        self.order = order
        number_of_elements = len(vector)
        dfdn = np.zeros(number_of_elements)
        
        if self.order == 1:
            
            for i in range(number_of_elements-1):
                if i  == 0:
                    f_forward = vector[i+1]
                    f_back = 0
                elif i == number_of_elements:
                    f_forward = 0
                    f_back = vector[i-1]
                else:
                    f_back = vector[i-1]
                    f_forward = vector[i+1]
                
                dfdn[i] = (f_forward - f_back) / dn
            
            return dfdn
        
        elif self.order == 2:
            
            for i in range(number_of_elements-1):
                if i  == 0:
                    f_forward = vector[i+1]
                    f_back = 0
                elif i == number_of_elements:
                    f_forward = 0
                    f_back = vector[i-1]
                else:
                    f_back = vector[i-1]
                    f_forward = vector[i+1]
                f = vector[i]
                dfdn[i] = (f_forward - 2 * f + f_back) / (dn ** 2)
            #print(dfdn)
            return dfdn
            
        else:
            raise ValueError('this order derivative is not implemented')
           
    def ForwardEuler(self, dc):
        
        #t_forward = t_current + dt * f_of_t_current
        pass
    
    def onedimension(self, c, dt, dx, L = 1, alpha = 1, beta = 1, Nt = 50):
        c = np.array(c)
        c_hold = []
        if self.method == 'forward':
            for i in range(Nt):
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * self.forward(c, dx, order = 2)            
                dc = -L * dFdc            
                c = c + dt * dc #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold
            return c
        
        elif self.method == 'backward':
            for i in range(Nt):
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * self.backward(c, dx, order = 2)            
                dc = -L * dFdc            
                c = c + dt * dc #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold
            return c
        
        elif self.method == 'central':   
            for i in range(Nt):
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * self.central(c, dx, order = 2)            
                dc = -L * dFdc            
                c = c + dt * dc    #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold
            return c_hold
 
    
    def twodimension(self, c, dt, dx, dy, L = 1, alpha = 1, beta = 1, Nt = 50):
        c = np.array(c)
        c_hold = []
        if self.method == 'forward':
            for i in range(Nt):
                grad_c = ( self.forward(c, dx, order = 2)+ self.forward(c, dy, order = 2))            
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * grad_c            
                dc = -L * dFdc            
                c = c + dt * dc #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold
            return c
        
        elif self.method == 'backward':
            for i in range(Nt):
                grad_c = ( self.backward(c, dx, order = 2)+ self.backward(c, dy, order = 2))            
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * grad_c            
                dc = -L * dFdc            
                c = c + dt * dc #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold
            return c
        
        elif self.method == 'central':
            for i in range(Nt):
                grad_c = ( self.central(c, dx, order = 2)+ self.central(c, dy, order = 2))         
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * grad_c          
                dc = -L * dFdc           
                c = c + dt * dc  #Forward Euler Method
            self.c_hold = c_hold
            return c
        
        
    
    def threedimension(self, c, dt, dx, dy, dz, L = 1, alpha = 1, beta = 1, Nt = 50):
        c = np.array(c)
        c_hold = []
        if self.method == 'forward':
            for i in range(Nt):
                grad_c = ( self.forward(c, dx, order = 2)+ self.forward(c, dy, order = 2) + self.forward(c, dz, order = 2))
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * grad_c            
                dc = -L * dFdc            
                c = c + dt * dc #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold
            return c
        
        elif self.method == 'backward':
            for i in range(Nt):
                grad_c = ( self.backward(c, dx, order = 2)+ self.backward(c, dy, order = 2) + self.backward(c, dz, order = 2))
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * grad_c            
                dc = -L * dFdc            
                c = c + dt * dc #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold
            return c
        
        elif self.method == 'central':
            for i in range(Nt):
                grad_c = ( self.central(c, dx, order = 2)+ self.central(c, dy, order = 2) + self.central(c, dz, order = 2))
                dFdc = 2 * alpha * c * (2 * c ** 2 - 3 * c + 1) - 2 * beta * grad_c            
                dc = -L * dFdc            
                c = c + dt * dc #Forward Euler Method
                c_hold.append(c)
            self.c_hold = c_hold            
            return c

'''
SquaredMash: Use to create a mesh grid with constant variations in x,y,z.

        Attribuites:
            number_of_elements: number of grid points in each direction. Default = 50
            
            
        Methods: 
            onedimensional: Creates a one dimension vector/array
            twodimensional: Creates a two dimension vector/array
            threedimensional: Creates a three dimension vector/array
    

'''

class SquaredMesh:
    def __init__(self, number_of_elements = 50):
        self.number_of_elements = number_of_elements
        
    
    def onedimension(self, start, stop, num = None):
        if not num:
            num = self.number_of_elements

        self.dx = abs((start - stop)/num)
        
        return np.linspace(start = start, stop = stop, num = num), self.dx
    
    def twodimension(self, x_start, x_stop, y_start = None, y_stop = None, num = None):
        if not num:
            num = self.number_of_elements
        if not y_start:
            y_start = x_start
            y_stop = x_stop
        
        self.dx = abs((x_start - x_stop)/num)
        self.dy = abs((y_start - y_stop)/num)
        
        x = np.linspace(start = x_start, stop = x_stop, num = num)
        y = np.linspace(start = y_start, stop = y_stop, num = num)
        
        mesh1, mesh2 = np.meshgrid(x,y)
        
        return mesh1*mesh2, self.dx, self.dy
    
    def threedimension(self, x_start, x_stop, y_start = None, y_stop = None,z_start = None, z_stop = None, num = None):
        if not num:
            num = self.number_of_elements
        if not y_start:
            y_start = x_start
            y_stop = x_stop
        if not z_start:
            z_start = x_start
            z_stop = x_stop
        
        self.dx = abs((x_start - x_stop)/num)
        self.dy = abs((y_start - y_stop)/num)
        self.dz = abs((z_start - z_stop)/num)
        
        x = np.linspace(start = x_start, stop = x_stop, num = num)
        y = np.linspace(start = y_start, stop = y_stop, num = num)
        z = np.linspace(start = z_start, stop = z_stop, num = num)
        
        mesh1, mesh2, mesh3 = np.meshgrid(x,y,z)
        
        return mesh1, self.dx, self.dy, self.dz
    
    
    
    
    
    
    
    
    