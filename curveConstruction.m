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
  v     = @(t) sqrt(a^2 + b^2);
  kappa = @(t) abs(a) / (a^2 + b^2);
  tau   = @(t) b / (a^2 + b^2);
  
  range = linspace(0,6*pi,1000);

  simpleHelix = simpleHelix3D(a,b)(range);
  
  e0  = [0, 1/sqrt(2), 1/sqrt(2)];
  n0  = [-1, 0, 0];
  b0  = [0, -1/sqrt(2), 1/sqrt(2)];
  fi0 = [1,0,0];
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

##### Calculations for simple helix #####

# fi = (a*cos(t), a*sin(t), b*t)
# fi' = (-a*sin(t), a*cos(t), b)
# fi'(0) = (0, a, b)
# |fi'(0)| = sqrt(a^2 + b^2)
# 
# a = b = 1
# e(0) = (0, 1/sqrt(2), 1/sqrt(2))
# 
# fi'' = (-a*cos(t), -a*sin(t), 0)
# fi''(0) = (-a, 0, 0)
# <fi''(0), e(0)> = 0*(-a) + 0*(1/sqrt(2)) + 0*(1/sqrt(2)) = 0
# n*(0) = (-1, 0, 0)
# n(0) = (-1, 0, 0) 
# 
# b(0) = e(0) * n(0) = 0*e1 - 1/sqrt(2) * e2  - 1/sqrt(2) * e3 = (0, -1/sqrt(2), 1/sqrt(2))
# 
# e1       e2       e3 
# 0    1/sqrt(2)  1/sqrt(2)
# -1       0         0 
# 
# e(0) = (0, 1/sqrt(2), 1/sqrt(2))
# n(0) = (-1, 0, 0) 
# b(0) = (0, -1/sqrt(2), 1/sqrt(2))