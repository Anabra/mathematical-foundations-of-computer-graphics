# Util
1;

# f : R -> [R;R]
function DS = diffCurve2(f, xs, n = 1)
  dXS = diff(xs);
  dYS = diff(f(xs),1,2);
  DS1 = dYS ./ dXS;
  
  if (n == 1)
    DS = DS1;
  elseif (n == 2)
    dDS1 = diff(DS1,1,2);
    DS2  = dDS1 ./ dXS(1:end-1);
    DS   = DS2;
  else 
    error("Cannot differentiate curve more than twice");
  endif
  
endfunction


function VS = v(f, xs)
  VS = ones(1,5);
  DS = diffCurve2(f, xs);
  for i = 1:length(DS)
    VS(i) = norm(DS(:,i));
  endfor
endfunction


function ES = e2(f, xs)
  VS = v(f,xs);
  DS = diffCurve2(f, xs);
  ES = DS ./ VS;
endfunction


function NS = n2(f, xs)
  ES = e2(f,xs);
  NS(1,:) = -1 * ES(2,:);
  NS(2,:) = ES(1,:); 
endfunction  


function KS = kappa2(f, xs) 
  ddf = diffCurve2(f, xs, 2);
  NS  = n2(f,xs);
  VS2 = v(f,xs) .^ 2;
  
  NS  = NS(:,1:end-1);
  VS2 = VS2(:,1:end-1);
  
  KS = dot(ddf, NS) ./ VS2;
endfunction

function VNS = vnorm(vs)
  VNS = sqrt(dot(vs,vs));
end

function [CPS, RS, tmp] = tangentCircle2(f, xs)
  FS  = f(xs);         # n
  KS  = kappa2(f, xs); # n-2
  ROS = 1 ./ KS;       # n-2
  
  df  = diffCurve2(f, xs, 1); # n-1
  df  = df(:,1:end-1);        # n-2
  ddf = diffCurve2(f, xs, 2); # n-2
  CS  = dot(ddf, df) ./ dot(df,df);
  
  ddf_s  = ddf .- (CS .* df);
  ddf_s0 = ddf_s ./ vnorm(ddf_s);
  # NS = n2(f, xs); # ddphi_s_0 for 2-dim case, n-1
  
  FS = FS(:,1:end-2);
  # NS = NS(:,end-1);
  
  CPS = FS .+ (ROS .* ddf_s0);
  RS  = ROS;
  tmp = ddf_s0;
endfunction


function plotCurve2(xsys)
  xs = xsys(1,:);
  ys = xsys(2,:);
  plot(xs,ys);
endfunction


function plotCircle2(x,y,r)
  ps = 0:0.01:(2*pi);
  xs = x + r*cos(ps);
  ys = y + r*sin(ps);
  plot(xs,ys); 
end


function plotTangentCircle2(f, xs)
  FS = f(xs);
  [CPS, RS, _] = tangentCircle2(f, xs);
  
  xMin = min([FS(1,:), CPS(1,:)]);
  xMax = max([FS(1,:), CPS(1,:)]);
  yMin = min([FS(2,:), CPS(2,:)]);
  yMax = max([FS(2,:), CPS(2,:)]);
  
  gblMin = min([xMin, yMin]) - 0.5;
  gblMax = max([xMax, yMax]) + 0.5;
  
  figure;
  axis([gblMin, gblMax, gblMin, gblMax], 'square');
  
  hold on;
  curvePlot = plot(FS(1,:),FS(2,:));  
  
  # comet(CPS(1,:), CPS(2,:), 0.001);
  
  for t = 1:length(CPS)
    hold off;
    plot(FS(1,:),FS(2,:));
    axis([gblMin, gblMax, gblMin, gblMax], 'square');
    hold on; 
    plot(CPS(1,1:t), CPS(2,1:t)); 
    plotCircle2(CPS(1,t),CPS(2,t),RS(t));
    pause(0.00005); 
  end
  
  hold off;
  
endfunction

########################### 3-dimensional ###########################

# f : R -> [R;R]
# TODO: implement this using recursion
function DS = diffCurve3(f, ts, n = 1)
  dTS = diff(ts); 
  DS = diffCurve3Helper(dTS, f(ts), n);  
end

# ts: parameter values 
# FS: 3D function values
function DS = diffCurve3Helper(dTS, FS, n)  
  if (n == 0)
    DS = FS;
  elseif (n > 0)
    dFS = diff(FS, 1, 2);
    
    df  = dFS ./ dTS;
    DS  = diffCurve3Helper(dTS(1:end-1), df, n-1);
  else 
    error("Cannot differentiate curve negative number of times");
  endif
  
end

function VS = v3(f, xs)
  DS = diffCurve3(f, xs);
  
  k  = length(DS);
  VS = ones(1,k);

  for i = 1:k
    VS(i) = norm(DS(:,i));
  end
end

# CPS: center points of the circles 
# RS: radii of the circles 
# ES: unit length derivatives
function [CPS, RS, ES, NS] = tangentCircle3(f, ts)
  dFS  = diffCurve3(f, ts, 1)(:, 1:end-1);
  ddFS = diffCurve3(f, ts, 2);
  
  ES = dFS ./ norm(dFS);
  NS = ddFS .- dot(ddFS, ES).*ES;
  
  VS  = v3(f, ts)(1:end-1);
  VS2 = VS .^ 2;
  KS  = dot(ddFS, NS) ./ VS2;
  ROS = 1 ./ KS;
  
  CS     = dot(ddFS, dFS) ./ dot(dFS,dFS);
  ddf_s  = ddFS .- (CS .* dFS);
  ddf_s0 = ddf_s ./ vnorm(ddf_s);
  
  FS = f(ts)(:,1:end-2);
  
  CPS = FS .+ (ROS .* ddf_s0);
  RS  = ROS;
end

function plotCurve3(curve)
  if (rows(curve) != 3)
    curve = curve';
  end
  
  xs = curve(1,:);
  ys = curve(2,:);
  zs = curve(3,:);
  
  plot3(xs,ys,zs);
  xlabel('x');
  ylabel('y');
  zlabel('z');
end

# cp, v1, v2 are 3D COLUMN vectors 
# v1 /= c*v2
# r is scalar
function plotCircle3(cp, r, v1, v2)
  ps = linspace(0, 2*pi, 200);
  ONB = orth([v1,v2]);
  V1 = ONB(:,1);
  V2 = ONB(:,2);
  circle = cp + r*(cos(ps) .* V1) + r*(sin(ps) .* V2);
  
  plot3(circle(1,:), circle(2,:), circle(3,:));
  xlabel('x');
  ylabel('y');
  zlabel('z');
end

function plotTangentCircle3(f, ts)
  FS = f(ts);
  [CPS, RS, ES, NS] = tangentCircle3(f, ts);
  
  # xMin = min([FS(1,:), CPS(1,:)]);
  # xMax = max([FS(1,:), CPS(1,:)]);
  # yMin = min([FS(2,:), CPS(2,:)]);
  # yMax = max([FS(2,:), CPS(2,:)]);
  
  # gblMin = min([xMin, yMin]) - 0.5;
  # gblMax = max([xMax, yMax]) + 0.5;
  
  figure;
  # axis([gblMin, gblMax, gblMin, gblMax], 'square');
  
  hold on;
  curvePlot = plot3(FS(1,:), FS(2,:), FS(3,:));  
  
  # comet(CPS(1,:), CPS(2,:), 0.001);
  
  for t = 1:length(CPS)
    hold off;
    plot3(FS(1,:), FS(2,:), FS(3,:));
    # axis([gblMin, gblMax, gblMin, gblMax], 'square');
    hold on; 
    plot3(CPS(1,1:t), CPS(2,1:t), CPS(3,1:t)); 
    plotCircle3(CPS(:,t), RS(t), ES(:,t), NS(:,t));
    pause(0.00005); 
  end
  
  hold off;
  
end

### Examples

# logSpiral = @(x) [(e.^x).*cos(x); (e.^x).*sin(x)];
function FUN = logSpiral3D(a,b)
  r   = @(ts) a*(e.^(b*ts));
  FUN = @(ts) [r(ts).*cos(ts); r(ts).*sin(ts); ts];
end

function FUN = helix3D(a)
  r   = @(ts) a*ts;
  FUN = @(ts) [r(ts).*cos(ts); r(ts).*sin(ts); ts];
end

# a, b are given constants
# v     = @(t) sqrt(a^2 + b^2)
# kappa = @(t) abs(a) / (a^2 + b^2)
# tau   = @(t) b / (a^2 + b^2) 
function FUN = simpleHelix3D(a,b)
  FUN = @(ts) [a*cos(ts); a*sin(ts); b*ts];
end

# v     = @(t) sqrt(1 + cos(t)^2)
# kappa = @(t) sqrt(5 + 3*cos(t)^2) / (1 + cos(t)^2)^(3/2)
# tau   = @(t) 6*cos(t) / (5 + 3*cos(t))
function FUN = viviani()
  FUN = @(ts) [cos(ts).^2; cos(ts).*sin(ts); sin(ts)];
end