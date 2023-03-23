function [x_opt,J,history,searchdir,exitflag,output] = runfmincon(simparams)

optoptions = simparams.optoptions;
x = simparams.x0;
P_initial = simparams.P_initial;
mu = simparams.mu;

% Set up shared variables with outfun
% returnhere = 1;
history.x = [];
history.J = [];
history.t = [];
searchdir = [];
history.stepsize = [];
history.gradient = [];
history.firstorderopt = [];
% history.trustregionradius = [];
history.constrviolation = [];
% t_last = x(end-1:end,:);

% Call optimization
optoptions.OutputFcn = @outfun;

[x_opt,J,exitflag,output] = fmincon(@(x)obj_min_tcm(x,simparams),x,[],[],[],[],[],[],@(x)constraint_min_tcm(x,simparams),optoptions)


     function stop = outfun(x,optimValues,state)
         stop = false;

         switch state
             case 'init'
                 hold on
             case 'iter'
             % Concatenate current point and objective function
             % value with history. x must be a row vector.
               history.J = [history.J; optimValues.fval];
               history.t = [history.t; toc];
               history.stepsize = [history.stepsize; optimValues.stepsize];
               history.gradient = [history.gradient, optimValues.gradient];
               history.firstorderopt = [history.firstorderopt; optimValues.firstorderopt];
%                history.trustregionradius = [history.trustregionradius; optimValues.trustregionradius];
               history.constrviolation = [history.constrviolation; optimValues.constrviolation];


               if length(history.x) == 0
                   history.x = x;
               else
                   history.x(:,:,end+1) = x;
               end

               if mod(size(history.x,3),10)==0
                   optimValues.fval
               end
             % Concatenate current search direction with 
             % searchdir.
%                searchdir = [searchdir;... 
%                             optimValues.searchdirection'];
%                plot(x(1),x(2),'o');
             % Label points with iteration number and add title.
             % Add .15 to x(1) to separate label from plotted 'o'.
%                text(x(1)+.15,x(2),... 
%                     num2str(optimValues.iteration));
%                title('Sequence of Points Computed by fmincon');
             case 'done'
                 hold off
             otherwise
         end
     end



end

