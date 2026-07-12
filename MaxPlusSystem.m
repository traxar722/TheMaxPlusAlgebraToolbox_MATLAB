classdef MaxPlusSystem
    properties
        A   % cell array of (n x n) maxplus matrices (state->state)
        B   % cell array of (n x m) maxplus matrices (input->state)
        C   % cell array of (p x n) maxplus matrices (state->output)
        D   % cell array of (p x m) maxplus matrices (input->output)



        % The following parameters exist only if the constructor's
        % parameter is a PetriNetwork object
        input_transitions  % column vector of transitions classified as input transitions
        output_transitions % column vector of transitions classified as output transitions
        state_transitions  % column vector of transitions classified as state transitions

        J   % positions state->state
        Js  % positions input->state
        Jo  % positions state->output
        Jos % positions input->output

        %maximums; % this vector stores the following:
        maxJ   % maximum number of tokens for positions state->state
        maxJs  % maximum number of tokens for positions input->state
        maxJo  % maximum number of tokens for positions state->output
        maxJos % maximum number of tokens for positions input->output

        A_plus; % Aplus matrix from the  PetriNetwork object
        A_minus;  % Aneg matrix from the  PetriNetwork object
        Marking;    % Initial marking from the PetriNetwork object
    end

    methods
        function obj = MaxPlusSystem(pn)
            if nargin == 0 || isempty(pn)


                obj.A = [];
                obj.B = [];
                obj.C = [];
                obj.D = [];
                return ;


            end

            if ~(isobject(pn) && isa(pn, 'Petrinetwork') && strcmp(pn.Mode,'ptimed'))
                error('The parameter entered is not a PetriNetwork object of ''ptimed'' type.');
            end

            if check_marked_graph(pn) ~= 1
                error('The PetriNetwork object entered is not a marked graph.')
            end

            Aplus = pn.Aplus;   % matrix Aplus from object
            Aminus = pn.Aminus; % matrix Aminus from object
            M0 = pn.marking(:); % initial marking from object
            obj.A_plus = Aplus;
            obj.A_minus = Aminus;
            obj.Marking = M0;


            Temporization = [];
            if isprop(pn, 'timedVec') && ~isempty(pn.timedVec)
                Temporization = pn.timedVec(:);
            end

            % transitions
            [obj.input_transitions, obj.output_transitions, obj.state_transitions] = pn.classify_transitions();

            % positions
            [obj.J, obj.Js, obj.Jo, obj.Jos] = pn.classify_positions();
            % max values of tokens
            [obj.maxJ, obj.maxJs, obj.maxJo, obj.maxJos] = pn.maximum_tokens();

            % sizes for the matrices of the max-plus system
            n = numel(obj.state_transitions);
            m_inputs = numel(obj.input_transitions);
            p_outputs = numel(obj.output_transitions);

            % limits representing the number of matrices of each type
            % (starting from 0 -> limit)
            limA = 0; if ~isempty(obj.maxJ),  limA = obj.maxJ;  end
            limB = 0; if ~isempty(obj.maxJs), limB = obj.maxJs; end
            limC = 0; if ~isempty(obj.maxJo), limC = obj.maxJo; end
            limD = 0; if ~isempty(obj.maxJos),limD = obj.maxJos; end

            % initializing cell arrays with zeros(...,'maxplus')
            obj.A = cell(1, limA + 1);
            for k = 0:limA
                obj.A{k+1} = zeros(n, n, 'maxplus'); % A = repmat(-Inf, n, n);
            end

            obj.B = cell(1, limB + 1);
            for k = 0:limB
                obj.B{k+1} = zeros(n, max(1, m_inputs), 'maxplus'); % B = repmat(-Inf, n, max(1, m_inputs));
            end
            if m_inputs == 0, obj.B = {}; end

            obj.C = cell(1, limC + 1);
            for k = 0:limC
                obj.C{k+1} = zeros(max(1, p_outputs), n, 'maxplus'); %C = repmat(-Inf,max(1, p_outputs), n);
            end
            if p_outputs == 0, obj.C = {}; end

            obj.D = cell(1, limD + 1);
            for k = 0:limD
                obj.D{k+1} = zeros(max(1, p_outputs), max(1, m_inputs), 'maxplus'); %D = repmat(-Inf,max(1, p_outputs),max(1, m_inputs));
            end
            if p_outputs == 0 || m_inputs == 0, obj.D = {}; end

            % getTiming function that returns Temporization(j) if valid,
            % otherwise 0 - NESTED function!!!
            function timing = getTiming(j)
                if ~isempty(Temporization) && numel(Temporization) >= j && j >= 1
                    timing = Temporization(j);
                else
                    timing = 0;
                end
            end

            % populate A from J (state->state)
            for idx = 1:numel(obj.J)
                place = obj.J(idx);
                k = M0(place);
                val = getTiming(place);

                t_src = find(Aplus(:, place) == 1, 1); % source
                t_dst = find(Aminus(:, place)== 1, 1); % destination

                col = find(obj.state_transitions == t_src, 1); %column
                lin = find(obj.state_transitions == t_dst, 1); %row

                % replace the zero from the maxplus matrix (zero=-inf) with
                % the corresponding timing value from the temporization
                % vector
                if ~isempty(lin) && ~isempty(col)
                    idxCell = k + 1;
                    if idxCell > numel(obj.A)
                        obj.A{idxCell} = zeros(n, n, 'maxplus');
                    end
                    obj.A{idxCell}(lin, col) = obj.A{idxCell}(lin, col) + maxplus(val);
                end
            end

            % populate B from Js (input->state)
            for idx = 1:numel(obj.Js)
                place = obj.Js(idx);
                k = M0(place);
                val = getTiming(place);

                t_src = find(Aplus(:, place) > 0, 1); %source
                t_dst = find(Aminus(:, place) > 0, 1); %destination

                col = find(obj.input_transitions == t_src, 1);  %column
                lin = find(obj.state_transitions == t_dst, 1);  %row

                % replace the zero from the maxplus matrix (zero=-inf) with
                % the corresponding timing value from the temporization
                % vector
                if ~isempty(lin) && ~isempty(col)
                    idxCell = k + 1;
                    if isempty(obj.B) || idxCell > numel(obj.B)
                        if m_inputs > 0, obj.B{idxCell} = zeros(n, m_inputs, 'maxplus'); end
                    end
                    obj.B{idxCell}(lin, col) = obj.B{idxCell}(lin, col) + maxplus(val);
                end
            end

            % populate C from Jo (state->output)
            for idx = 1:numel(obj.Jo)
                place = obj.Jo(idx);
                k = M0(place);
                val = getTiming(place);

                t_src = find(Aplus(:, place) > 0, 1); %source
                t_dst = find(Aminus(:, place) > 0, 1); %destination

                col = find(obj.state_transitions == t_src, 1); %column
                lin = find(obj.output_transitions == t_dst, 1); %row

                %replace the zero from the maxplus matrix (zero=-inf) with
                %the corresponding timing value from the temporization
                %vector
                if ~isempty(lin) && ~isempty(col)
                    idxCell = k + 1;
                    if isempty(obj.C) || idxCell > numel(obj.C)
                        if p_outputs > 0, obj.C{idxCell} = zeros(p_outputs, n, 'maxplus'); end
                    end
                    obj.C{idxCell}(lin, col) = obj.C{idxCell}(lin, col) + maxplus(val);
                end
            end

            % populate D from Jos (input->output)
            for idx = 1:numel(obj.Jos)
                place = obj.Jos(idx);
                k = M0(place);
                val = getTiming(place);

                t_src = find(Aplus(:, place) > 0, 1); %source
                t_dst = find(Aminus(:, place) > 0, 1); %destination

                col = find(obj.input_transitions == t_src, 1); %column
                lin = find(obj.output_transitions == t_dst, 1); %row

                % replace the zero from the maxplus matrix (zero=-inf) with
                % the corresponding timing value from the temporization
                % vector
                if ~isempty(lin) && ~isempty(col)
                    idxCell = k + 1;
                    if isempty(obj.D) || idxCell > numel(obj.D)
                        if p_outputs > 0 && m_inputs > 0
                            obj.D{idxCell} = zeros(p_outputs, m_inputs, 'maxplus');
                        end
                    end
                    obj.D{idxCell}(lin, col) = obj.D{idxCell}(lin, col) + maxplus(val);
                end
            end
        end

        function [X_mat, Y_mat] = sim(mpsystem, steps_or_u)
            % sim - Simulates a MaxPlusSystem object over N steps.
            %
            %   [X, Y] = sim(mpsystem, N)          % autonomous system (N = number of steps)
            %   [X, Y] = sim(mpsystem, U)          % non-autonomous system (U = m x N maxplus input matrix)
            %
            %   X_mat : n x N  matrix, each column is x(k)
            %   Y_mat : p x N  matrix, each column is y(k)  (empty if there are no outputs)

            if ~(isobject(mpsystem) && isa(mpsystem, 'MaxPlusSystem'))
                error('The first parameter must be a MaxPlusSystem object');
            end

            % number of each tyoe of transition
            n = numel(mpsystem.state_transitions);
            m = numel(mpsystem.input_transitions);
            p = numel(mpsystem.output_transitions);

            is_autonomous = isempty(mpsystem.B) && isempty(mpsystem.D); % or there are no input transitions


            if is_autonomous
                if ~isscalar(steps_or_u) || mod(steps_or_u, 1) ~= 0 || steps_or_u < 0
                    error('For autonomous systems, provide a nonnegative integer number of steps.');
                end
                N = steps_or_u;
                U = [];
            else
                % U must be provided as a matrix (m x N), where m is the number
                % of input transitions and each column represents
                % u(1)...u(k) input vectors of the system
                if isscalar(steps_or_u)
                    error('The second parameter must be a matrix (m x N) where each column is the k-th firing time of each input transition');
                end
                U = steps_or_u;

                %make the matrix into a maxplus object
                if ~isa(U, 'maxplus')
                    U = maxplus(U);
                end

                %verify if the U matrix is of m x N dimensions where m is
                %the number of  input transitions and N is the number of
                %iterations wanted
                if size(U, 1) ~= m
                    error('U input matrix must have %d rows (one for each input transition)', m);
                end
                % N is the number of iterations and is equal to the number
                % of columns of U (matrix)
                N = size(U, 2);
            end

            %verify if the MaxPlusSystem object has atleast A0 matrix
            %corresponding to the state transitions
            if ~isempty(mpsystem.A) && numel(mpsystem.A) >= 1
                if (is_nill(mpsystem.A{1})~=1)
                    error('The system entered has a DEADLOCK');
                end
                A0_star = kstar(mpsystem.A{1}); % using the kstar function implemented in the maxplus class
            else
                A0_star = eye(n, 'maxplus');   %the kstar matrix is equal to the identity matrix in maxplus algebra
            end


            X_mat = zeros(n,N,'maxplus');   % state  matrix
            Y_mat = zeros(max(1,p),N,'maxplus');  % output matrix


            for k = 1:N


                % for each state transition t_i, x_bar_i(k) is equal to 0 if the transition is validated for the k-th time(if every
                % input place of t_i has at least k tokens in the initial
                % marking M0)
                % otherwise x_bar_i(k) = epsilon (-inf).
                % x_bar must have be of dimensions n x 1 where n is the
                % number of state transitions
                x_bar = zeros(n,1,'maxplus') ;  % initialise to epsilon (-inf)

                for i = 1:n
                    t_idx = mpsystem.state_transitions(i);


                    % returns the indexes of the lines that have elements >
                    % 0 ( positions that go into the transition and are needed to validate it)
                    input_places = find(mpsystem.A_minus(t_idx, :) > 0);

                    if isempty(input_places) %this means that the respective transition is a input transition
                        % input transitions are always validated
                        x_bar(i) = maxplus(0);
                    else
                        % enabled at step k only if every input place has >= k tokens
                        if all(mpsystem.Marking(input_places) >= k) %the transition is validated for the k-th iteration if there are atleast k tokens in the initial marking
                            x_bar(i) = maxplus(0);
                        end
                    end
                end

                %calculating v_k
                v_k = x_bar;   %intialized with x_bar(k)


                for l = 1:(numel(mpsystem.A) - 1)
                    %A matrixes
                    Al     = mpsystem.A{l + 1};
                    x_prev = get_history(X_mat, k - l, n);
                    v_k    = v_k + Al * x_prev;
                end

                % B matrixes
                if ~is_autonomous && ~isempty(mpsystem.B)
                    for ls = 0:(numel(mpsystem.B) - 1)
                        Bls   = mpsystem.B{ls + 1};
                        u_cur = get_input(U, k - ls, m);
                        v_k   = v_k + Bls * u_cur;
                    end
                end

                % calculate the state equation x(k) = A0_star * v(k)
                X_mat(:, k) = A0_star * v_k;


                % y_bar(k) is the analogue of x_bar for output transitions: it is epsilon when
                % there are no output transitions, so we initialise to epsilon
                if p > 0

                    y_bar = zeros(p,1,'maxplus');   % initialise to epsilon (-inf)

                    for i = 1:p
                        t_idx = mpsystem.output_transitions(i);

                        input_places = find(mpsystem.A_minus(t_idx, :) > 0);

                        if isempty(input_places)
                            y_bar(i) = maxplus(0);
                        else
                            if all(mpsystem.Marking(input_places) >= k)
                                y_bar(i) = maxplus(0);
                            end
                        end
                    end

                    y_k = y_bar;   %initialiazed with y_bar(k)

                    % C matrixes
                    if ~isempty(mpsystem.C)
                        for lr = 0:(numel(mpsystem.C) - 1)
                            Clr    = mpsystem.C{lr + 1};
                            x_prev = get_history(X_mat, k - lr, n);
                            y_k    = y_k + Clr * x_prev;
                        end
                    end

                    % D matrixes
                    if ~is_autonomous && ~isempty(mpsystem.D)
                        for lrs = 0:(numel(mpsystem.D) - 1)
                            Dlrs  = mpsystem.D{lrs + 1};
                            u_cur = get_input(U, k - lrs, m);
                            y_k   = y_k + Dlrs * u_cur;
                        end
                    end

                    Y_mat(:, k) = y_k;
                end
            end

            % if there are  no outputs, return empty
            if p == 0
                Y_mat = [];
            end



            function x = get_history(X_matrix, idx, num_states)
                % returns x(idx); for idx <= 0 returns the epsilon vector (-inf).
                if idx <= 0
                    x = zeros(num_states,1,'maxplus');   % epsilon = -inf
                else
                    x = X_matrix(:, idx); %column
                end
            end

            function u = get_input(U_mat, idx, num_inputs)
                % returns u(idx); for idx <= 0 or empty returns the epsilon vector.
                if idx <= 0 || isempty(U_mat)
                    u = maxplus(repmat(-Inf, num_inputs, 1));   % epsilon = -inf
                else
                    u = U_mat(:, idx); %column
                end
            end
            function plot_sim(X, Y, U)
                % steps k=1..N (N from X). Each row uses a distinct color;

                if nargin < 1, error('At least X must be provided'); end
                if nargin < 2 || isempty(Y), Y = []; end
                if nargin < 3 || isempty(U), U = []; end

                N = size(X,2);        % reference number of time steps for X
                tX = 1:N;

                if ~isempty(Y)
                    if size(Y,2) ~= N
                        error('Y must have the same number of columns (time steps) as X.');
                    end
                    tY = tX;
                end
                if ~isempty(U)
                    if size(U,2) ~= N
                        error('U must have the same number of columns (time steps) as X.');
                    end
                    tU = tX;
                end

                % markers and sizes
                mX = '*'; szX = 60;
                mY = 'v'; szY = 60;
                mU = 'x'; szU = 60;

                % plot X: one color per row
                nX = size(X,1);
                colsX = lines(max(1,nX));
                figure('Name','Plot for X'); hold on; grid on;
                for i = 1:nX
                    c = colsX(mod(i-1,size(colsX,1))+1, :);
                    scatter(tX, X(i,:), szX, 'Marker', mX, 'MarkerEdgeColor', c);
                end
                xlim([0, N + 1]); xticks(0:1:N+1);
                xlabel('k'); ylabel('value'); title('States (X)');
                legend("i=" + (1:nX), 'Location', 'bestoutside');
                hold off;

                % plot Y
                if ~isempty(Y)
                    nY = size(Y,1);
                    colsY = lines(max(1,nY));
                    figure('Name','Plot for Y'); hold on; grid on;
                    for i = 1:nY
                        c = colsY(mod(i-1,size(colsY,1))+1, :);
                        scatter(tY, Y(i,:), szY, 'Marker', mY, 'MarkerEdgeColor', c);
                    end
                    xlim([0, N + 1]); xticks(0:1:N+1);
                    xlabel('k'); ylabel('value'); title('Outputs (Y)');
                    legend("i=" + (1:nY), 'Location', 'bestoutside');
                    hold off;
                end

                % plot U
                if ~isempty(U)
                    nU = size(U,1);
                    colsU = lines(max(1,nU));
                    figure('Name','Plot for U'); hold on; grid on;
                    for i = 1:nU
                        c = colsU(mod(i-1,size(colsU,1))+1, :);
                        scatter(tU, U(i,:), szU, 'Marker', mU, 'MarkerEdgeColor', c);
                    end
                    xlim([0, N + 1]); xticks(0:1:N+1);
                    xlabel('k'); ylabel('value'); title('Inputs (U)');
                    legend("i=" + (1:nU), 'Location', 'bestoutside');
                    hold off;
                end
            end



            X_double = double(X_mat);
            Y_double = double(Y_mat);
            U = double(U);
            plot_sim(X_double,Y_double,U);
        end

     function [inTransPerMachine, outTransPerMachine] = getMachineTransitions(obj, machinePlaces)
            % getMachineTransitions  Get places to their input/output transitions per machine.
            % returns local row indices into the X matrix produced by sim()
            % roth indices are looked up first in state_transitions (rows of X), then in
            % output_transitions (rows of Y), and returned as negative values for output
            % transitions so they can be differentiated.
            %
            % Inputs:
            %   machinePlaces - (MxK) matrix of nonneg integer place indices (0 used as padding to keep matrix dimensions).
            % Outputs:
            %   inTransPerMachine  - (MxK) local row index in X (or Y if negative) for the
            %                        transition that puts tokens into the place.
            %   outTransPerMachine - (MxK) local row index in X (or Y if negative) for the
            %                        transition that takes tokens from the place.

            if ~(isobject(obj) && isa(obj, 'MaxPlusSystem'))
                error('This method must be called on a MaxPlusSystem object.');
            end
            if isempty(obj.A_plus) || isempty(obj.A_minus)
                error('Object must contain A_plus and A_minus matrices.');
            end

            Aplus  = obj.A_plus;
            Aminus = obj.A_minus;

            validateattributes(machinePlaces, {'numeric'}, {'2d', 'nonempty'});
            if any(mod(machinePlaces(:),1) ~= 0) || any(machinePlaces(:) < 0)
                error('machinePlaces must contain nonnegative integers (zeros allowed as padding).');
            end

            [M, K] = size(machinePlaces);
            [~, P] = size(Aplus);

            nonzeroIdx = machinePlaces(:) > 0;
            if any(machinePlaces(nonzeroIdx) > P)
                error('machinePlaces contains place indices larger than number of places (%d).', P);
            end

            inTransPerMachine  = zeros(M, K);
            outTransPerMachine = zeros(M, K);

            for mi = 1:M
                for kj = 1:K
                    place = machinePlaces(mi, kj);
                    if place == 0
                        continue;
                    end

                    % global transition index that generates tokens INTO this place (A_plus col)
                    global_in  = find(Aplus(:, place)  > 0, 1);
                    % global transition index that consumes tokens FROM this place (A_minus col)
                    global_out = find(Aminus(:, place) > 0, 1);

                    % mapping global_in to local row in X (state_transitions) or Y (output_transitions)
                    % positive value = row in X, negative value = row in Y
                    if ~isempty(global_in)
                        local_in = find(obj.state_transitions == global_in, 1);
                        if isempty(local_in)
                            local_in_y = find(obj.output_transitions == global_in, 1);
                            if ~isempty(local_in_y)
                                local_in = -local_in_y; % mark as Y-row
                            else
                                local_in = 0;
                            end
                        end
                        inTransPerMachine(mi, kj) = local_in;
                    end

                    % mapping global_out to local row in X (state_transitions) or Y (output_transitions)
                    if ~isempty(global_out)
                        local_out = find(obj.state_transitions == global_out, 1);
                        if isempty(local_out)
                            local_out_y = find(obj.output_transitions == global_out, 1);
                            if ~isempty(local_out_y)
                                local_out = -local_out_y;
                            else
                                local_out = 0;
                            end
                        end
                        outTransPerMachine(mi, kj) = local_out;
                    end
                end
            end
        end


        function plotMachineGantt(obj, machinePlaces, X, Y)
            % plotMachineGantt  Draw per-machine Gantt charts from transition firing times.
            %
            % Inputs:
            %   obj           - MaxPlusSystem object
            %   machinePlaces - (MxK) matrix of nonneg integer place indices (0 used as padding)
            %   X             - (n x N) matrix of state-transition firing times from sim()
            %   Y             - (p x N) matrix of output-transition firing times from sim()
            %                   (optional if there are no output
            %                   transitions)
            %
            % For each place j and each step k
            % the interval that the token spends from being deposited to being consumed is:
            %   start  = firing time of the transition that deposits the token into j at step k
            %   finish = firing time of the transition that consumes the token from j at step k
            % both start and finish refer to the SAME step k.

            barWidth = 0.75;

            % input validation
            if ~(isobject(obj) && isa(obj, 'MaxPlusSystem'))
                error('Must call on a MaxPlusSystem object.');
            end
            if isempty(obj.A_plus) || isempty(obj.A_minus)
                error('Object must contain A_plus and A_minus matrices.');
            end
            validateattributes(machinePlaces, {'numeric'}, {'2d','nonempty'});
            if any(mod(machinePlaces(:),1) ~= 0) || any(machinePlaces(:) < 0)
                error('machinePlaces must contain nonnegative integers (zeros allowed).');
            end

            % X must have n rows (state transitions) and not T (all transitions)
            n = numel(obj.state_transitions);
            if size(X, 1) ~= n
                error('X must have %d rows (one per state transition).', n);
            end
            N = size(X, 2);
            if N < 1
                error('X must have at least one column.');
            end

            % Y must have p rows and not T
            if nargin < 4 || isempty(Y)
                Y = zeros(numel(obj.output_transitions), N);
            end
            p = numel(obj.output_transitions);
            if size(Y, 1) ~= p
                error('Y must have %d rows (one per output transition).', p);
            end

            %  get firing time for a local index (positive = X row, negative = Y row)
            function t = getFiringTime(localIdx, k)
                if localIdx > 0
                    t = X(localIdx, k);
                elseif localIdx < 0
                    t = Y(-localIdx, k);
                else
                    t = -Inf;
                end
            end

            % get local transition indices
            [inTransPerMachine, outTransPerMachine] = obj.getMachineTransitions(machinePlaces);

            M    = size(machinePlaces, 1);
            cmap = lines(M);

            figure;
            for mi = 1:M
                subplot(M, 1, mi);
                hold on;

                mp    = machinePlaces(mi, :);
                rects = zeros(0, 3); % [k, startT, finishT]

                for kj = 1:numel(mp)
                    place  = mp(kj);
                    if place == 0, continue; end

                    inIdx  = inTransPerMachine(mi, kj);
                    outIdx = outTransPerMachine(mi, kj);
                    if inIdx == 0 || outIdx == 0, continue; end

                    for k = 1:N
                        % start = firing of the input transition at step k
                        % finish = firing of the output transition at the SAME step k
                        startT  = getFiringTime(inIdx,  k);
                        finishT = getFiringTime(outIdx, k);

                        if ~isfinite(startT) || ~isfinite(finishT)
                            continue;
                        end
                        if finishT <= startT
                            continue;
                        end
                        rects(end+1, :) = [k, startT, finishT];
                    end
                end

                if isempty(rects)
                    title(sprintf('Machine %d: no tasks to plot', mi));
                    ylim([0 1]);
                    hold off;
                    continue;
                end

                % one bar per row in rects, y position = row index
                nTasks = size(rects, 1);
                for t = 1:nTasks
                    startT  = rects(t, 2);
                    finishT = rects(t, 3);
                    yBot = t - barWidth/2;
                    yTop = t + barWidth/2;

                    % correct vertex order for a patch rectangle: DL to DR to TR to TL
                    px = [startT; finishT; finishT; startT];
                    py = [yBot;   yBot;    yTop;    yTop  ];
                    patch(px, py, cmap(mi,:), 'EdgeColor', 'k', 'FaceAlpha', 0.85);

                    % label k for each bar : step number centered inside the bar
                    text((startT+finishT)/2, t, sprintf('k=%d', rects(t,1)), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment',   'middle', ...
                        'FontSize', 8, 'Color', 'k');
                end

                ylim([0.5, nTasks + 0.5]);
                yticks(1:nTasks);
                xlabel('Time');
                ylabel('Task');
                title(sprintf('Resource %d', mi));
                grid on;
                hold off;
            end
        end





    end
end