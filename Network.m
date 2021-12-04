%Load all variables


clear variables;


prompt = '\nPlease provide name of input file for a transportation problem: ';
x = input(prompt, 's');


start = tic;
fid = fopen(x);

if fid < 0
    fprintf('\nInput file cannot be opened...exiting\n');
    return
end

oname = strcat('Output', x);
fileID = fopen(oname,'w');

tline = fgetl(fid); %Skip description line
tline = fgetl(fid);
%Extract the dimension from TopLine
n = str2double(tline);

tline = fgetl(fid); %Skip description line
%populate b
tline = fgetl(fid);
names = split(tline);
b = cellfun(@str2double,names);

tline = fgetl(fid); %Skip description line

C = zeros(n, n);

%populate C
for i=1:n
   tline = fgetl(fid);
   namei = split(tline);
   constrainti = cellfun(@str2double,namei);
   for j=1:n
       C(i, j) = constrainti(j);
   end
end

fclose(fid); %Close Input

fprintf(fileID,'The Network problem is:\n');
transposeb = b.';

fprintf(fileID, "\nb = \n\t");
for ii = 1:size(transposeb,1)
    fprintf(fileID,'%g\t',transposeb(ii,:));
    fprintf(fileID,'\n\t');
end

fprintf(fileID, "\nC = \n\t");
for ii = 1:size(C,1)
    fprintf(fileID,'%g\t',C(ii,:));
    fprintf(fileID,'\n\t');
end
clear tranposeb;

%Check if balanced
if sum(b) ~= 0
    fprintf('\nProblem is unbalanced\n');
    fprintf(fileID, '\nProblem is unbalanced\n');
    return
end




%Add the Central Warehouse Node
n1 = n+1;
%Create a structure to show edges (n+1 x n+1 matrix)
A = zeros(n1,n1);

%Set the initial edges to Central Warehouse
BasicVariables = zeros(n1, n1);
%Blocked is -1
for i=1:n
    if b(i) > 0
        BasicVariables(i, n1) = 1;
        BasicVariables(n1, i) = -1;
        A(i, n1) = b(i);
    else
        BasicVariables(i, n1) = -1;
        BasicVariables(n1, i) = 1;
        A(n1, i) = -b(i);
    end
end
BasicVariables(n1, n1) = -1;

fprintf('\nStarting Phase I\n');
fprintf(fileID,'\nStarting Phase I\n');

iteration = 0;

c = BasicVariables;

%Phase I
while 1
    iteration = iteration + 1;

    fprintf('\nSimplex Iteration %d\n', iteration);
    fprintf(fileID, '\nSimplex Iteration %d\n', iteration);
    
    %Obtain lambdai and lambdaj
    lambdai = zeros(n1, 1);
    lambdaj = zeros(n1, 1);
    
    iFilled = zeros(n1,1);
    jFilled = zeros(n1,1);
    iFilled(n1) = 1;
    jFilled(n1) = 1;

    while 1
        
        nochange = 0;

        bt = 0;
        for i=1:n1
            for j=1:n1
                if BasicVariables(i,j) == 1
                    if c(i,j) ~= lambdai(i) - lambdaj(j)
                        bt = 1;
                        break;
                    end
                end
            end
            if bt == 1
                break
            end
        end
        
        if bt == 0
            break
        end
        

        for i=1:n1
            for j=1:n1

                if BasicVariables(n1-i+1,n1-j+1) == 1
                    if iFilled(n1-i+1) == 0 && jFilled(n1-j+1) == 1
                        lambdai(n1-i+1) = c(n1-i+1, n1-j+1) + lambdaj(n1-j+1);
                        iFilled(n1-i+1) = 1;
                        nochange = 1;
                    end
                    if jFilled(n1-j+1) == 0 && iFilled(n1-i+1) == 1
                        lambdaj(n1-j+1) =  lambdai(n1-i+1) - c(n1-i+1, n1-j+1);
                        jFilled(n1-j+1) = 1;
                        nochange = 1;
                    end
                end
            end
        end 

        if nochange == 0
            for i=1:n1
                if iFilled(i) == 0 && jFilled(i) == 1
                    lambdai(i) = lambdaj(i);
                    iFilled(i) = 1;
                elseif jFilled(i) == 0 && iFilled(i) == 1
                    lambdaj(i) = lambdai(i);
                    jFilled(i) = 1;
                elseif iFilled(i) && jFilled(i)
                    continue
                end
            end
        end
    end

   
    
    for i=1:n1
        if iFilled(i) == 0 && jFilled(i) == 1
            lambdai(i) = lambdaj(i);
            iFilled(i) = 1;
        elseif jFilled(i) == 0 && iFilled(i) == 1
            lambdaj(i) = lambdai(i);
            jFilled(i) = 1;
        elseif iFilled(i) == 0 && jFilled(i) == 0
            lambdai(i) = -1;
            lambdaj(i) = lambdai(i);
        elseif iFilled(i) && jFilled(i)
            continue
        end
    end
    
    r = zeros(n1,n1);
    
    %Fill in r for non-basic
    for i=1:n1
        for j = 1:n1
            if BasicVariables(i,j) == 0
                r(i,j) = c(i,j) - (lambdai(i) - lambdaj(j));
            end
        end
    end
    
    %Find an rij<0
   
    bt = 0;
    smallest = 0;
    for i=1:n1
        for j=1:n1
            if r(i,j) < 0
                iB = i;
                jB = j;
                smallest = r(i,j);
                bt = 1;
                break
            end
        end
        if bt ==1
            break
        end
    end

    if smallest>=0   %Optimal
        fprintf('\nInitial BFS Reached\n\n');
        fprintf(fileID, '\nInitial BFS Reached\n\n');
        break
    end
    
    fprintf('\nX%d%d will enter basis\n',iB,jB);
    fprintf(fileID, '\nX%d%d will enter basis\n',iB,jB);
    BasicVariables(iB,jB) = 1;
    
    %Construct Graph
    s = [];
    t = [];
    for i=1:n1
        for j=1:n1
            if BasicVariables(i,j) == 1
                s = [s i];
                t = [t j];
            end
        end
    end
    G = graph(s,t);
    %Obtain Cycle
    cycle = allcycles(G);
    cycle = cell2mat(cycle);
    
    cyclelength = size(cycle, 2);
    %Check if second element is correct, if not then all the values will
    %have a multiplier of -1
    
    %Start from our iB
    for i=1:cyclelength
        if cycle(i) == iB
            break
        end
    end
    startingindex = i;

    newcycle = zeros(cyclelength,1);
    for i=1:cyclelength
        newcycle(i) = cycle(rem(((i-1) + (startingindex-1)), cyclelength) + 1);
    end
    %Now we have correct start
    cycle = [newcycle; iB];
    %Must make sure order is correct
    if cycle(2) ~= jB
        cycle = flip(cycle);
    end



    %get edges 
    s = [];
    t = [];

    for i=1:cyclelength
        s = [s cycle(i)];
        t = [t cycle(i+1)];
    end

    %determine negative edges and positive edges
    %there are cyclelength edges
    
    signvector = zeros(n1, n1);

    thetavector = [];
    for i=1:cyclelength
        if BasicVariables(t(i), s(i)) == 1
            signvector(t(i), s(i)) = -1; 
            thetavector = [thetavector A(t(i), s(i))];
        elseif BasicVariables(s(i), t(i)) == 1
            signvector(s(i), t(i)) = 1;
        else
            continue
        end
    end
    theta = min(thetavector);


    %Find the basic variable 
    for i=1:cyclelength
        if A(t(i), s(i)) == theta && signvector(t(i),s(i)) == -1 
            dI = t(i);
            dJ = s(i);
            break;
        end
    end
    
    %update along cycle
    for i=1:n1
        for j=1:n1
            if signvector(i,j) == -1
                A(i,j) = A(i,j) - theta;
            elseif signvector(i,j) == 1
                    A(i,j) = A(i,j) + theta;
            else
                continue
            end
        end
    end
    
    
    fprintf('X%d%d will leave basis\n',dI,dJ);
    fprintf(fileID, 'X%d%d will leave basis\n',dI,dJ);
    
    BasicVariables(dI, dJ) = 0;
    
end


fprintf("Initial BFS : \n");
fprintf(fileID, "\nInitial BFS = \n");
for i=1:n
    for j=1:n
        if BasicVariables(i,j) == 1
            fprintf('X%d%d = %f, ', i, j, A(i,j));
            fprintf(fileID, 'X%d%d = %f, ', i, j, A(i,j));
        end
    end
end
fprintf('\n');
fprintf(fileID, '\n');




BasicVariables = BasicVariables(1:n, 1:n);
A = A(1:n, 1:n);
c = C;
n1 = n;



z = sum(sum(A.*c));
fprintf('Initial Cost = %f', z);
fprintf(fileID, 'Initial Cost = %f', z);



fprintf('\n');
fprintf(fileID, '\n');


fprintf('\nStarting Phase II\n');
fprintf(fileID, '\nStarting Phase II\n');

iteration = 0;

%Phase II
while 1
    iteration = iteration + 1;

    
    fprintf('\nSimplex Iteration %d\n', iteration);
    fprintf(fileID, '\nSimplex Iteration %d\n', iteration);
    
    %Obtain lambdai and lambdaj
    lambdai = zeros(n1, 1);
    lambdaj = zeros(n1, 1);
    
    iFilled = zeros(n1,1);
    jFilled = zeros(n1,1);
    iFilled(n1) = 1;
    jFilled(n1) = 1;

    while 1
        
        nochange = 0;

        bt = 0;
        for i=1:n1
            for j=1:n1
                if BasicVariables(i,j) == 1
                    if c(i,j) ~= lambdai(i) - lambdaj(j)
                        bt = 1;
                        break;
                    end
                end
            end
            if bt == 1
                break
            end
        end
        
        if bt == 0
            break
        end
        

        for i=1:n1
            for j=1:n1

                if BasicVariables(n1-i+1,n1-j+1) == 1
                    if iFilled(n1-i+1) == 0 && jFilled(n1-j+1) == 1
                        lambdai(n1-i+1) = c(n1-i+1, n1-j+1) + lambdaj(n1-j+1);
                        iFilled(n1-i+1) = 1;
                        nochange = 1;
                    end
                    if jFilled(n1-j+1) == 0 && iFilled(n1-i+1) == 1
                        lambdaj(n1-j+1) =  lambdai(n1-i+1) - c(n1-i+1, n1-j+1);
                        jFilled(n1-j+1) = 1;
                        nochange = 1;
                    end
                end
            end
        end 

        if nochange == 0
            for i=1:n1
                if iFilled(i) == 0 && jFilled(i) == 1
                    lambdai(i) = lambdaj(i);
                    iFilled(i) = 1;
                elseif jFilled(i) == 0 && iFilled(i) == 1
                    lambdaj(i) = lambdai(i);
                    jFilled(i) = 1;
                elseif iFilled(i) && jFilled(i)
                    continue
                end
            end
        end
    end

   
    
    for i=1:n1
        if iFilled(i) == 0 && jFilled(i) == 1
            lambdai(i) = lambdaj(i);
            iFilled(i) = 1;
        elseif jFilled(i) == 0 && iFilled(i) == 1
            lambdaj(i) = lambdai(i);
            jFilled(i) = 1;
        elseif iFilled(i) == 0 && jFilled(i) == 0
            lambdai(i) = -1;
            lambdaj(i) = lambdai(i);
        elseif iFilled(i) && jFilled(i)
            continue
        end
    end
    
    r = zeros(n1,n1);
    
    %Fill in r for non-basic
    for i=1:n1
        for j = 1:n1
            if BasicVariables(i,j) == 0
                r(i,j) = c(i,j) - (lambdai(i) - lambdaj(j));
            end
        end
    end
    
    %Find an rij<0
   
    bt = 0;
    smallest = 0;
    for i=1:n1
        for j=1:n1
            if r(i,j) < 0
                iB = i;
                jB = j;
                smallest = r(i,j);
                bt = 1;
                break
            end
        end
        if bt ==1
            break
        end
    end

    if smallest>=0   %Optimal
        fprintf('\nOptimal Solution Reached\n\n');
        fprintf(fileID,'\nOptimal Solution Reached\n\n');
        break
    end
    
    fprintf('\nX%d%d will enter basis\n',iB,jB);
    fprintf(fileID, '\nX%d%d will enter basis\n',iB,jB);
    BasicVariables(iB,jB) = 1;
    
    %Construct Graph
    s = [];
    t = [];
    for i=1:n1
        for j=1:n1
            if BasicVariables(i,j) == 1
                s = [s i];
                t = [t j];
            end
        end
    end
    G = graph(s,t);
    %Obtain Cycle
    cycle = allcycles(G);
    cycle = cell2mat(cycle);
    
    cyclelength = size(cycle, 2);
    %Check if second element is correct, if not then all the values will
    %have a multiplier of -1
    
    %Start from our iB
    for i=1:cyclelength
        if cycle(i) == iB
            break
        end
    end
    startingindex = i;

    newcycle = zeros(cyclelength,1);
    for i=1:cyclelength
        newcycle(i) = cycle(rem(((i-1) + (startingindex-1)), cyclelength) + 1);
    end
    %Now we have correct start
    cycle = [newcycle; iB];
    %Must make sure order is correct
    if cycle(2) ~= jB
        cycle = flip(cycle);
    end



    %get edges 
    s = [];
    t = [];

    for i=1:cyclelength
        s = [s cycle(i)];
        t = [t cycle(i+1)];
    end

    %determine negative edges and positive edges
    %there are cyclelength edges
    
    signvector = zeros(n1, n1);

    thetavector = [];
    for i=1:cyclelength
        if BasicVariables(t(i), s(i)) == 1
            signvector(t(i), s(i)) = -1; 
            thetavector = [thetavector A(t(i), s(i))];
        elseif BasicVariables(s(i), t(i)) == 1
            signvector(s(i), t(i)) = 1;
        else
            continue
        end
    end
    theta = min(thetavector);


    %Find the basic variable 
    for i=1:cyclelength
        if A(t(i), s(i)) == theta && signvector(t(i),s(i)) == -1 
            dI = t(i);
            dJ = s(i);
            break;
        end
    end
    
    %update along cycle
    for i=1:n1
        for j=1:n1
            if signvector(i,j) == -1
                A(i,j) = A(i,j) - theta;
            elseif signvector(i,j) == 1
                    A(i,j) = A(i,j) + theta;
            else
                continue
            end
        end
    end
    
    fprintf('X%d%d will leave basis\n',dI,dJ);
    fprintf(fileID, 'X%d%d will leave basis\n',dI,dJ);
    BasicVariables(dI, dJ) = 0;

    fprintf("\nCurrent BFS : \n");
    fprintf(fileID, "\nCurrent BFS = \n");
        for i=1:n
            for j=1:n
                if BasicVariables(i,j) == 1
                    fprintf('X%d%d = %f, ', i, j, A(i,j));
                    fprintf(fileID, 'X%d%d = %f, ', i, j, A(i,j));
                end
            end
        end

        fprintf('\n');
        fprintf(fileID, '\n');

        z = sum(sum(A.*c));
        fprintf('Current Cost = %f\n', z);
        fprintf(fileID, 'Current Cost = %f\n', z);
    
end

fprintf("Optimal Solution : \n");
fprintf(fileID, "\nOptimal Solution = \n");
for i=1:n
    for j=1:n
        if BasicVariables(i,j) == 1
            fprintf('X%d%d = %f, ', i, j, A(i,j));
            fprintf(fileID, 'X%d%d = %f, ', i, j, A(i,j));
        end
    end
end



fprintf('\n');
fprintf(fileID, '\n');

z = sum(sum(A.*c));
fprintf('Optimal Cost = %f\n', z);
fprintf(fileID, 'Optimal Cost = %f\n', z);



fprintf('\nContraint Check:\n');
fprintf(fileID, '\nContraint Check:\n');

allpos = 0;
for i=1:n
    for j=1:n
        if A(i,j) < 0
            fprintf('X%d%d < 0, violation\n', i, j);
            fprintf(fileID, 'X%d%d < 0, violation\n', i, j);
            allpos = 1;
        end
    end
end

if allpos == 0
    fprintf('All Xij >= 0\n');
    fprintf(fileID,'All Xij >= 0\n');
end

anyvio = 0;

%Checking as sources
fprintf('\nChecking Sources\n');
fprintf(fileID, '\nChecking Sources\n');
for i=1:n
    vio = 0;
    if b(i) > 0
        fprintf('Source %d, b(%d) and all recieved= %f, sum(A(%d)) = %f, ',i, i, b(i) + sum(A(:,i)), i, sum(A(i, :)));
        fprintf(fileID, 'Source %d, b(%d) and all recieved= %f, sum(A(%d)) = %f, ',i, i, b(i) + sum(A(:,i)), i, sum(A(i, :)));
        if  b(i) + sum(A(:,i)) ~= sum(A(i, :))
            vio = 1;
        end

        if vio == 1
            fprintf('Violation.\n');
            fprintf(fileID, 'Violation.\n');
            anyvio =1;
        else
            fprintf('Correct.\n');
            fprintf(fileID, 'Correct.\n');
        end
    end
end


fprintf('\nChecking Destination\n');
fprintf(fileID, '\nChecking Destination\n');
for i=1:n
    vio = 0;
    if b(i) <= 0
        fprintf('Destination %d, b(%d) and all lost= %f, sum(A(%d)) = %f, ',i, i, -b(i) + sum(A(i,:)), i, sum(A(:, i)));
        fprintf(fileID, 'Destination %d, b(%d) and all lost= %f, sum(A(%d)) = %f, ',i, i, -b(i) + sum(A(i,:)), i, sum(A(:, i)));
        if  b(i) + sum(A(:,i)) ~= sum(A(i, :))
            vio = 1;
        end

        if vio == 1
            fprintf('Violation.\n');
            fprintf(fileID, 'Violation.\n');
            anyvio =1;
        else
            fprintf('Correct.\n');
            fprintf(fileID, 'Correct.\n');
        end
    end

end


if anyvio == 1 || allpos == 1
    fprintf('\nThe solution is infeasible.');
    fprintf(fileID, '\nThe solution is infeasible.');
else
    fprintf('\nThe solution is feasible.');
    fprintf(fileID, '\nThe solution is feasible.');
end

Elapsed_time = toc(start);
fprintf(fileID, "\n\nTime Taken: %f seconds\n", Elapsed_time);
fprintf("\n\nTime Taken: %f seconds\n", Elapsed_time);
fclose(fileID);

clear variables;

