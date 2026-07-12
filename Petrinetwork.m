classdef Petrinetwork
    properties
        Aplus           % matrix (nonnegative integer elements)
        Aminus          % matrix (nonnegative integer elements)
        marking         % vector with length == number of columns of A
        timedVec        % optional vector dependent of the selected mode (size depends on Mode)
        Mode = ''       % '' | 'untimed' | 'ptimed' | 'ttimed'
    end


    methods
        function obj = Petrinetwork(Aplus, Aminus, marking, varargin)
            % Petrinetwork(Aplus, Aminus, marking, 'Mode', mode, 'TimedVec', timedVec)
            if nargin < 3
                error('Provide Aplus, Aminus and marking.');
            end

            validateattributes(Aplus, {'numeric'}, {'2d'});
            validateattributes(Aminus, {'numeric'}, {'2d'});
            if size(Aplus) ~= size(Aminus)
                error('Aplus and Aminus must have the same dimensions.');
            end

            % all elements of the matrixes need to be >= 0 and integer
            if any(mod(Aplus(:),1) ~= 0) || any(mod(Aminus(:),1) ~= 0) || any(Aplus(:) < 0) || any(Aminus(:) < 0)
                error('All elements of Aplus and Aminus must be nonngegative integers.');
            end

            [m,n] = size(Aplus);
            %size of the marking vector must match the number of columns of
            %the matrixes
            validateattributes(marking, {'numeric'}, {'vector', 'numel', n});
            marking = marking(:); % column

            % optional name-value pairs
            p = inputParser;
            addParameter(p, 'Mode', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'TimedVec', [], @(x) isnumeric(x) || isempty(x));
            parse(p, varargin{:});
            mode = char(lower(p.Results.Mode));
            timedVec = p.Results.TimedVec;


            % validate mode and timedVec requirements


            if isempty(mode) || strcmp(mode, 'untimed')
                % untimed and empty mode do not require timedVec
                timedVec = [];
                if isempty(mode)
                    mode = '';
                end
            elseif strcmp(mode, 'ptimed')
                % timedVec required and size must match the number of
                % column of the matrixes
                if isempty(timedVec)
                    error('TimedVec is required for mode ''ptimed''.');
                end
                validateattributes(timedVec, {'numeric'}, {'vector', 'numel', n});
                timedVec = timedVec(:);
            elseif strcmp(mode, 'ttimed')
                % timedVec required and size must match the number of
                % rows of the matrixes
                if isempty(timedVec)
                    error('TimedVec is required for mode ''ttimed''.');
                end
                validateattributes(timedVec, {'numeric'}, {'vector', 'numel', m});
                timedVec = timedVec(:);
            else
                error('Mode must be '''' or ''untimed'' or ''ptimed'' or ''ttimed''.');
            end

            % assign
            obj.Aplus = Aplus;
            obj.Aminus = Aminus;
            obj.marking = marking;
            obj.Mode = mode;
            if ~isempty(timedVec)
                obj.timedVec = timedVec;
            end
        end

        function verify = check_marked_graph(petrinet)
            if ~(isobject(petrinet) && isa(petrinet, 'Petrinetwork'))
                error('The parameter entered is not a PetriNetwork object.');
            end

            A_plus = petrinet.Aplus;
            A_minus = petrinet.Aminus;
            % verify that the matrixes Aplus and Aminus have all elements 0
            % or 1
            if max(max(abs(A_plus - A_minus))) >1
                verify=0;
                return
            end

            [~,n]=size(A_plus);
            % verify that the number of input transitions is equal to the
            % number of output transitions from a place ( this number needs
            % to be exactly 1 for each type of transition)
            for j = 1:n
                if sum(A_plus(:,j)) >1  % sum(A_plus(:,j)) = number of input transitions in place j
                    verify=0;
                    return
                end
                if sum(A_minus(:,j)) >1   % sum(A_mins(:,j)) = number of output transitions from place j
                    verify=0;
                    return
                end
            end
            verify=1; %returns 1 if the petrinet is a marked graph
            return
        end


        function [input_transitions, output_transitions, state_transitions] = classify_transitions(petrinet)
            % Classify transitions of a Petrinetwork object
            %   - Returns row indices of transitions that are:
            %       input_transitions  : have the value 1 in Aplus
            %       matrix("positive" arcs)
            %       output_transitions : have the value 1 in Aminus matrix
            %       ("negative" arcs)
            %       state_transitions  : have the value one in both Aplus
            %       and Aminus matrixes

            if ~(isobject(petrinet) && isa(petrinet, 'Petrinetwork') && strcmp(petrinet.Mode,'ptimed'))
                error('The parameter entered is not a PetriNetwork object of ''ptimed'' type.');
            end

            if ~(check_marked_graph(petrinet))
                error('The Petri Network entered is not a marked graph.');
            end

            A_plus = petrinet.Aplus;
            A_minus = petrinet.Aminus;



            % detect presence of positive/negative arcs per row (transition)
            hasPos = any(A_plus > 0, 2);   %  row has at least one positive arc
            hasNeg = any(A_minus > 0, 2);   % row has at least one "negative" arc (Aminus>0)

            input_transitions  = find(hasPos & ~hasNeg);   % positives only
            output_transitions = find(hasNeg & ~hasPos);   % negatives only
            state_transitions  = find(hasPos & hasNeg);    % both
        end

        function [J, Js, Jo, Jos] = classify_positions(petrinet)

            %Classify positions of a Petrinetwork object
            %- Returns column indices of positions situated between the
            %following :
            %   J    : state -> state
            %   Js   : input -> state
            %   Jo   : state -> output
            %   Jos  : input -> output

            A_plus=petrinet.Aplus;
            A_minus= petrinet.Aminus;
            [~, n] = size(A_plus);

            J   = []; Js  = []; Jo  = []; Jos = [];
            [input,output,state]=classify_transitions(petrinet);
            for j = 1:n
                pos_rows = find(A_plus(:,j) > 0);
                neg_rows = find(A_minus(:,j) > 0);

                if isempty(pos_rows) || isempty(neg_rows)
                    continue %skip to the end of the for loop if there the vectors are empty
                end


                pos_idx = pos_rows(1);
                neg_idx = neg_rows(1);

                isPosInput  = any(pos_idx == input);
                isPosState  = any(pos_idx == state);
                isNegState  = any(neg_idx == state);
                isNegOutput = any(neg_idx == output);

                if isPosState && isNegState
                    J(end+1,1) = j; %adds an element to this vector initially blank
                elseif isPosInput && isNegState
                    Js(end+1,1) = j; %adds an element to this vector initially blank
                elseif isPosState && isNegOutput
                    Jo(end+1,1) = j; %adds an element to this vector initially blank
                elseif isPosInput && isNegOutput
                    Jos(end+1,1) = j; %adds an element to this vector initially blank
                end
            end
        end

        function [maxJ, maxJs, maxJo, maxJos] = maximum_tokens(petrinet)

            % Returns the maximum number of tokens for each type of position; empty if group has no indices.
            initial_marking = petrinet.marking;
            if isempty(initial_marking)
                error('marcaj must be a nonempty vector');
            end

            %column vector for indexing
            initial_marking = initial_marking(:);
            [J, Js, Jo, Jos] = classify_positions(petrinet);
            maxJ   = []; maxJs  = [];   maxJo= []; maxJos = [];

            if ~isempty(J)
                maxJ = max(initial_marking(J));
            end
            if ~isempty(Js)
                maxJs = max(initial_marking(Js));
            end
            if ~isempty(Jo)
                maxJo = max(initial_marking(Jo));
            end
            if ~isempty(Jos)
                maxJos = max(initial_marking(Jos));
            end
        end






    end
end
