# Constructing biregular 3D curves from v (velocity), kappa (curveture), tau (torsion)
2;

function fi = constructCurve(v, kappa, tau, range, e0 = [1,0,0], n0 = [0,1,0], b0 = [0,0,1], fi0 = [0,0,0]) 
  ODESFun = @(t,y) mkDiffEqSystem(t, y, v, kappa, tau);
  [_,y]   = ode45(ODESFun, range, [e0, n0, b0, fi0]);
  fi      = y(:, 10:12);
end

function ODEs = mkDiffEqSystem(t, y, v, kappa, tau)
  ODEs(1) = v(t)*kappa(t)*y(4);
  ODEs(2) = v(t)*kappa(t)*y(5);
  ODEs(3) = v(t)*kappa(t)*y(6);
  
  ODEs(4) = -v(t)*kappa(t)*y(1) + v(t)*tau(t)*y(7);
  ODEs(5) = -v(t)*kappa(t)*y(2) + v(t)*tau(t)*y(8);
  ODEs(6) = -v(t)*kappa(t)*y(3) + v(t)*tau(t)*y(9);
  
  ODEs(7) = -v(t)*tau(t)*y(4);
  ODEs(8) = -v(t)*tau(t)*y(5);
  ODEs(9) = -v(t)*tau(t)*y(6);
  
  ODEs(10) = v(t)*y(1);
  ODEs(11) = v(t)*y(2);
  ODEs(12) = v(t)*y(3);
end

function testSimpleHelix(a,b)
  c     = a^2 + b^2;
  
  v     = @(t) sqrt(c);
  kappa = @(t) abs(a) / (c);
  tau   = @(t) b / (c);
  
  range = linspace(0,6*pi,1000);

  simpleHelix = simpleHelix3D(a,b)(range);
  
  e0  = [0, a/sqrt(c), b/sqrt(c)];
  n0  = [-1, 0, 0];
  b0  = [0, -b/sqrt(c), a/sqrt(c)];
  fi0 = [a,0,0];
  curve = constructCurve(v, kappa, tau, range, e0, n0, b0, fi0)';
  
  # simpleHelix = simpleHelix + (curve(:,1) - simpleHelix(:,1));
  # curve       = rotateZ(pi/2) * rotateY(3*pi/4) * curve;
  
  #1 curve       = curve .+ [1;0;3*pi];
  #1 curve       = rotateZ(-pi) * rotateY(pi) * curve;
  #2 curve       = rotate(0, -pi/4, pi/2) * curve;
  
  plotCurve3(curve);
  hold on;
  plotCurve3(simpleHelix);
end

function testViviani()
  v     = @(t) sqrt(1 + cos(t)^2);
  kappa = @(t) sqrt(5 + 3 * cos(t)^2) / (1 + cos(t)^2)^(3/2);
  tau   = @(t) (6*cos(t)) / (5 + 3 * cos(t)^2);
  
  range = linspace(0,2*pi,1000);

  vivianiCurve = viviani()(range);
  
  e0  = [0, 1/sqrt(2), 1/sqrt(2)];
  n0  = [-1, 0, 0];
  b0  = [0, -1/sqrt(2), 1/sqrt(2)];
  fi0 = [1,0,0];
  curve = constructCurve(v, kappa, tau, range, e0, n0, b0, fi0)';
  
  plotCurve3(curve);
  hold on;
  plotCurve3(vivianiCurve);
end

function R = rotateX(t)
  R = [ 1,      0,       0
      ; 0, cos(t), -sin(t)
      ; 0, sin(t),  cos(t) 
      ];
end 

function R = rotateY(t)
  R = [ cos(t),  0, sin(t)
      ;      0,  1,      0
      ; -sin(t), 1, cos(t) 
      ];
end 

function R = rotateZ(t)
  R = [ cos(t), -sin(t), 0
      ; sin(t),  cos(t), 0
      ;      0,       0, 1 
      ];
end 

function R = rotate(t1, t2, t3) 
  R = rotateZ(t3) * rotateY(t2) * rotateX(t1);
end

################################## examples ####################################

# testSimpleHelix(1,1);

##### Calculations for simple helix #####

# fi = (a*cos(t), a*sin(t), b*t)
# fi' = (-a*sin(t), a*cos(t), b)
# fi'(0) = (0, a, b)
# |fi'(0)| = sqrt(a^2 + b^2)
# 
# c = a^2 + b^2
# a = b = 1
# e(0) = (0, a/sqrt(c), b/sqrt(c))
# 
# fi'' = (-a*cos(t), -a*sin(t), 0)
# fi''(0) = (-a, 0, 0)
# <fi''(0), e(0)> = 0*(-a) + 0*(a/sqrt(c)) + 0*(b/sqrt(c)) = 0
# n*(0) = (-1, 0, 0)
# n(0) = (-1, 0, 0) 
# 
# b(0) = e(0) * n(0) = 0*e1 - b/sqrt(c) * e2  + a/sqrt(c) * e3 = (0, -b/sqrt(c), a/sqrt(c))
# 
# e1       e2       e3 
# 0    a/sqrt(c)  b/sqrt(c)
# -1       0         0 
# 
# e(0) = (0, a/sqrt(c), a/sqrt(c))
# n(0) = (-1, 0, 0) 
# b(0) = (0, -b/sqrt(2), a/sqrt(c))





##### Calculations for Viviani curve #####

# fi   = (cos(t)^2, cos(t)*sin(t), sin(ts))
# fi'  = (-2*sin(t)*cos(t), cos(t)^2 - sin(t)^2, cos(t))
# fi'' = (-2*cos(2t), -2*sin(2t), -sin(t))

# fi'(0)   = (0, 1, 1)
# |fi'(0)| = sqrt(2)

# e(0)    = (0, 1/sqrt(2), 1/sqrt(2))
# fi''(0) = (-2, 0, 0) 
# <fi''(0), e(0)> = (0,0,0)
# n*(0) = fi''(0) - <fi''(0), e(0)> * e(0) = (-2,0,0)
# n(0) = (-1,0,0)

# b(0) = 0*e1 - (1/sqrt(2))*e2 + (1/sqrt(2))*e3 

# e1       e2       e3 
# 0    1/sqrt(2)  1/sqrt(2)
# -1       0         0 

# e(0) = (0, 1/sqrt(2), 1/sqrt(2))
# n(0) = (-1,0,0)
# b(0) = (0, -1/sqrt(2), 1/sqrt(2))