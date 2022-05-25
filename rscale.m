         function[Nbar]=rscale(a,b,c,d,k)
         % Given the single-input linear system:
         %       .
         %       x = Ax + Bu
         %       y = Cx + Du
         % and the feedback matrix K,
         %
         % the function rscale(sys,K) or rscale(A,B,C,D,K)
         % finds the scale factor N which will
         % eliminate the steady-state error to a step reference
         % for a continuous-time, single-input system
         % with full-state feedback using the schematic below:
         %
         %                         /---------\
         %      R         +     u  | .       |
         %      ---> N --->() ---->| X=Ax+Bu |--> y=Cx ---> y
         %                -|       \---------/
         %                 |             |
         %                 |<---- K <----|
         %
         %8/21/96 Yanjie Sun of the University of Michigan
         %        under the supervision of Prof. D. Tilbury
         %6/12/98 John Yook, Dawn Tilbury revised
         narginchk(2,5);
         % --- Determine which syntax is being used ---
         nargin1 = nargin;
         if (nargin1==2)	% System form
         		[A,B,C,D] = ssdata(a);
         		K=b;
         elseif (nargin1==5) % A,B,C,D matrices
         		A=a; B=b; C=c; D=d; K=k;
         else 
             error('Input must be of the form (sys,K) or (A,B,C,D,K)')
         end
         % compute Nbar
         s1 = size(A);
         s2 = size(B);
         Z = [zeros(s1(1),s1(2)+s2(2)); eye(s1(1),s1(2)+s2(2))];
         N = [A-eye(s1),B;C,D]\Z;
         Nx = N(1:s1(1),:);
         Nu = N(s1:end,:);
         Nbar=Nu + K*Nx;