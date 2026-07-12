classdef maxplus < double

    methods

        function obj = maxplus(a)
            % maxplus class constructor.
            if ~isreal(a), error('Max-plus matrices must be real.'), end
            obj = obj@double(a);
        end

        function z = plus(x,y)
            % maxplus addition

            sx = size(double(x));
            sy = size(double(y));
            if ~(isequal(sx, sy) || isscalar(x) || isscalar(y))
                error('Operands must have the same size.');
            end
            z = maxplus(max(double(x), double(y)));
        end

        function z = times(x,y)
            % maxplus scalar multiplication
            z = maxplus(double(x) + double(y));
        end

        function z = mtimes(x,y)
            % maxplus matrix multiplication
            [m,n] = size(x);
            [n1,p] = size(y);
            if n ~= n1, error('Matrix dimensions not compatible.'), end
            if max([m n p]) == 1
                % scalar multiplication
                z = maxplus(double(x) + double(y));
            else
                % matrix multiplication
                z = zeros(m,p,'maxplus');
                for i = 1:m
                    for j = 1:p
                        z(i,j) = maxplus(max(double(x(i,:).') + double(y(:,j))));
                    end
                end
            end
        end

        function z = sum(x, dim)
            % SUM  that calculates the sum on rows/columns/all elements


            if nargin < 2
                dim = 1;
            end

            % scalar case
            if isscalar(x)
                z = x;
                return;
            end

            % "all" case (all elements)
            if ischar(dim) && strcmp(dim, "all")
                z = zeros(1,1,'maxplus'); % zero element in max-plus
                for i = 1:numel(x)
                    z = z + x(i);
                end
                return;
            end

            % one dimensional array (vector) case
            if isvector(x)
                z = zeros(1,1,'maxplus');
                for i = 1:numel(x)
                    z = z + x(i);
                end
                return;
            end

            % two dimensional array (matrix) case
            [m, n] = size(x);
            if dim == 1 % calculate the sum on columns
                z = zeros(1,n,'maxplus') ;
                for j = 1:n
                    for i = 1:m
                        z(j) = z(j) + x(i, j);
                    end
                end
            elseif dim == 2 % calculate the sum on rows
                z = zeros(m,1,'maxplus') ;
                for i = 1:m
                    for j = 1:n
                        z(i) = z(i) + x(i, j);
                    end
                end
            else
                error('Second argument must be 1, 2, or "all".');
            end
        end

        function z = mpower(x,k)
            % MPOWER - matrix or scalar power
            % for scalar x, one can use either
            %         x^k with k a nonnegative integer
            % or
            %         x^(1/k) with k a positive integer
            %
            % for x a square matrix, k must be a nonnegative integer

            % scalar case
            
if isscalar(x)
    
    if ~( (mod(k,1) == 0 && k >= 0) || (k > 0 && abs(round(1/k) - 1/k) < eps) )
        error('Exponent must be a nonnegative integer, or of the form 1/k with k a positive integer for scalar.');
    end

    if mod(k,1) == 0
        
        if k == 0
            z = maxplus(0);
        else
            z = maxplus(double(x) * double(k));
        end
        return;
    else
        
        p = round(1/k);
        z = maxplus(double(x) / double(p));
        return;
    end
end

            % matrix case
            [m,n] = size(x);
            if m ~=n, error('Matrix is not square'), end
            if mod(k,1) ~= 0  || k < 0, error('Exponent must be nonnegative for matrix.'), end

            z = eye(size(x), 'maxplus');

            for i = 1:k
                z = z * x;
            end
        end



        function z = trace(x)
            [m,n] = size(x);
            if nargin < 1
                error('There must be atleast one input parameter.');
            end

            if m ~=n, error('Matrix is not square'), end
            z = sum(diag(x));
        end



        function nill = is_nill(x)

            if nargin < 1
                error('There must be atleast one input parameter.');
            end

            [m,n] = size(x);
            if m ~=n, error('Matrix is not square'), end

            % Verify if A^m is equal to the null matrix in max plus algebra


            if ((x^m)==zeros(m,'maxplus'))
                nill = 1;
            else
                nill = 0;
            end
        end


        function kst = kstar(x)
             if (is_nill(x) == 0)
                error('Matrix is not nilpotent')
            end
            [m, ~] = size(x);

            kst = eye(m, 'maxplus');
            for i = 1:m-1
                kst = kst + (x^i);
            end
        end


        function z = diag(v, k)
            % DIAG creates a diagonal matrix or extracts a diagonal from a matrix in max-plus algebra.
            % DIAG(V,K) when V is a vector with N components creates a square matrix
            % of order N + ABS(K) with the elements of V on the K-th diagonal.
            % K = 0 is the main diagonal, K > 0 is above the main diagonal, and K < 0 is below the main diagonal.
            %
            % DIAG(V) is the same as DIAG(V,0) and puts V on the main diagonal.
            %
            % DIAG(X,K) when X is a matrix returns a column vector formed from
            % the elements of the K-th diagonal of X.
            %
            % DIAG(X) returns the main diagonal of X. DIAG(DIAG(X)) is a diagonal matrix.
            if nargin < 2
                k = 0; %default main diagonal
            end

            if isvector(v)
                n = length(v);
                z = zeros(n + abs(k), n + abs(k),'maxplus');

                v_obj = maxplus(v); %make v into a maxplus object

                if k >= 0
                    for i = 1:n
                        z(i, i + k) = v_obj(i);
                    end
                else
                    for i = 1:n
                        z(i - k, i) = v_obj(i);
                    end
                end
            else
                [m, n] = size(v);
                if k >= 0
                    if k > (n - 1), error('Cannot extract diagonal'), end
                    z = zeros(1, min(m, n - abs(k)), 'maxplus');
                    for i = 1:min(m, n - k)
                        z(i) = v(i, i + k);
                    end
                else
                    if abs(k) > (m - 1), error('Cannot extract diagonal'), end
                    z = zeros(1, min(m-abs(k), n), 'maxplus');
                    for i = 1:min(m + k, n)
                        z(i) = v(i - k, i);
                    end
                end
            end
        end

        function lambda = max_eig(A)

            % lambda = SUM ( trace(A^i)^(1/i) )
            if nargin < 1
                error('There must be atleast one input parameter.');
            end

            [m, n] = size(A);
            if m ~= n, error('Matrix must be square'); end

            lambda = zeros(1,1,'maxplus') ;

            for i = 1:n
                Ai = A^i;
                tr = trace(Ai);
                lambda = lambda + tr^(1/i);
            end
        end

        function verify = verify_acyclic_graph(A)

            %   return 1 if precedence graph associated with A is acyclic, 0 otherwise.
            %   if the precedence graph associated with A is acyclic it
            %   means that there are no eigenvalues
            %   return 0 otherwise
            %   input A must be a square matrix
            if nargin < 1
                error('There must be atleast one input parameter.');
            end

            if isempty(A)
                error('The input parameter must not be empty');
            end

            [m,n] = size(A);
            if m ~= n
                error('The input parameter must be a square matrix');
            end

            %check if there is an arc connecting a node to itself of length 1
            if trace(A) > zeros(1,1,'maxplus');
                verify = 0;
                return
            end

            %check if there is an arc connecting a node to itself of length > 2
            for i = 2:m
                if trace(A^i) > zeros(1,1,'maxplus');
                    verify = 0;
                    return
                end
            end
            verify = 1;

        end


        function verify = verify_connected_graph(A)
   %   return 1 if precedence graph G(A) is connected
   %   return 0 otherwise
   %   input A must be a square matrix

   if nargin < 1
      error('There must be atleast one input parameter.');
   end
   if isempty(A)
      error('The input parameter must not be empty');
   end
   [m, n] = size(A);
   if m ~= n
      error('The input parameter must be a square matrix');
   end

   % transform the directed graph in a undirected graph
   for i = 1:m
      A(i, i) = 0;
      for j = 1:m
         if A(j, i) > zeros(1,1,'maxplus')
            A(j, i) = 0;
            A(i, j) = 0;
         end
      end
   end

   % calculate the transitive closure matrix representing if there is a
   % possibility of reaching a node from another node with any number of
   % arcs
   A_closure = A;
   A_power   = A;
   for i = 2:m
      A_power   = A_power * A;
      A_closure = A_closure + A_power;
   end

   % check if the initial node can reach any other node
   for i = 1:m
      if A_closure(i, 1) == zeros(1,1,'maxplus')
         verify = 0;
         return
      end
   end
   verify = 1;
end

        function verify = verify_str_connected_graph(A)

            %     return 1 if precedence graph G(A) is strongly connected wich also means:
            %     matrix A is irreducible
            %     matrix A has a single eigenvalue
            %     return  0 otherwise.
            %     input A must be a square matrix.
            if nargin < 1
                error('There must be atleast one input parameter.');
            end

            if isempty(A)
                error('The input parameter must not be empty. ');
            end
            [m,n] = size(A);
            if m ~= n
                error('The input parameter must be a square matrix');
            end


            B = A;
            for i = 2:m
                B =   B + A^i;
            end

            if any(B(:)== zeros(1,1,'maxplus'))
                verify = 0;
                return
            end
            verify = 1;
        end

        function [x, err] = solve_Axb(A, b)


            %       returns the largest x that is restricted by  Ax <= b
            %       if A and B are scalars, x = A-b

            if nargin < 2
                error('There must be atleast two input parameter.');
            end

            if (isempty(A) || isempty(b))
                error('The input variables must not be empty')
            end

            Id = ones(size(b));
            x = Id - ((Id - b)'*  A)';
           
            %verify if the equation has a solution that satisfies Ax = b or if the resulted x is
            %a subsolution
            if A*x == b
                err = 0;
            else
                err = 1;
            end
        end



        function [x, err] = solve_xAxb(A, b)
    if nargin < 2
        error('There must be at least two input parameters.');
    end

    if (isempty(A) || isempty(b))
        error('The input variables must not be empty')
    end

    x = kstar(A) * b;

    % verify the solution
  if all(x == A*x + b)
        err = 0;  
    else
        err = 1; 
    end
end

        function disp(a)
            disp(double(a));
        end
    
    end
    methods (Static)
        function A = zeros(m,n)
            if nargin == 1
                if isscalar(m), n = m; else, n = m(2); m = m(1); end
            end

            if nargin == 0
                A = maxplus(-inf);
            else
                A = maxplus(repmat(-inf,[m,n]));
            end
        end

        function A = eye(m,n)
            if nargin == 1
                if isscalar(m), n = m; else, n = m(2); m = m(1); end
            end
            if nargin == 0
                A = maxplus(0);
            else
                A = zeros(m,n,'maxplus');
                for i = 1:m
                    for j = 1:n
                        if(i==j)
                            A(i,j) = 0;
                        end
                    end
                end
            end
        end
    end
end