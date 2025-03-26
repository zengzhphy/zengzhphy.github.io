clear;clc;

%Preparation
n=10000;                         % Markov chain length: 10000;  sample size: 50000
ns=20;                          % 20×20 lattice
T_mc = [0.50:0.05:2.00, 2.05:0.01:2.45, 2.50:0.05:5.00];  % Temperature from 1 to 5
tic;                            % Timer, usually runs for about 2 minutes when n=5000

%Different temperature
for jj=1:1:size(T_mc,2)
    X=sign(rand(ns,ns));        % Initialize all spins in the same direction +1 (simulating starting from 0 Kelvin)
    % Markov chain length: 10000 Monte Carlo step
    for j=1:1:n 
        % one Monte Carlo step is one update of all Ising fields (= ns*ns updates)
        for row = 1:ns          % Sequentially select each lattice site, row and column stored in index[1,2]
            for col = 1:ns
                index = [row, col];  
                % Compute the row and column indices of the four nearest-neighboring sites
                % using periodic boundary conditions (PBC)
                % +1 to transform 0(20)~19 to 1~20
                tmp1=rem(index(1),ns)+1; tmp2=rem(index(1)+1,ns)+1; tmp3=rem(index(1)-1,ns)+1;
                tmp4=rem(index(2),ns)+1; tmp5=rem(index(2)+1,ns)+1; tmp6=rem(index(2)-1,ns)+1;
                % Compute the possible energy change
                cen=X(tmp1,tmp4);
                right=X(tmp1,tmp5);  left=X(tmp1,tmp6);
                up=X(tmp2,tmp4);     down=X(tmp3,tmp4);
                deE=2*cen*(right+left+up+down);
                % Decide whether to flip the spin based on the Metropolis criterion
                if rand < exp(-deE./T_mc(jj))
                    X(tmp1,tmp4)=-X(tmp1,tmp4);
                end
            end
        end
    end    


    % Sampling phase: 50,000 samples, still applying the Metropolis criterion to maintain equilibrium
    for j=1:1:5*n
        for row = 1:ns          % Sequentially select each lattice site, row and column stored in index[1,2]
            for col = 1:ns
                index = [row, col];  
                % Compute the row and column indices of the four nearest-neighboring sites
                % using periodic boundary conditions (PBC)
                % +1 to transform 0(20)~19 to 1~20
                tmp1=rem(index(1),ns)+1; tmp2=rem(index(1)+1,ns)+1; tmp3=rem(index(1)-1,ns)+1;
                tmp4=rem(index(2),ns)+1; tmp5=rem(index(2)+1,ns)+1; tmp6=rem(index(2)-1,ns)+1;
                % Compute the possible energy change
                cen=X(tmp1,tmp4);
                right=X(tmp1,tmp5);  left=X(tmp1,tmp6);
                up=X(tmp2,tmp4);     down=X(tmp3,tmp4);
                deE=2*cen*(right+left+up+down);
                % Decide whether to flip the spin based on the Metropolis criterion
                if rand < exp(-deE./T_mc(jj))
                    X(tmp1,tmp4)=-X(tmp1,tmp4);
                end
            end
        end
        % Compute the average magnetization for a specific configuration
        m(j)=abs(mean(mean(X))); 
        % Compute the average energy for a specific configuration
        Xt1=X;Xt1(1,:)=[];Xt1=[Xt1; X(1,:)];
        Xt2=X;Xt2(:,1)=[];Xt2=[Xt2, X(:,1)];
        e(j)=-mean(mean(X.*Xt1+X.*Xt2));
    end

    % Compute statistical quantities at a specific temperature
    m_bar(jj)=mean(m);   
    e_bar(jj)=mean(e);
    cv_bar(jj)= ns^2 * std(e)^2 ./ T_mc(jj)^2;
end
toc;

% Plot results for visualization
% Plot 1: Average Magnetization vs Temperature
figure(1);
plot(T_mc,m_bar,'k.', 'MarkerSize', 10);
xlabel('Temperature T (J)');  
ylabel('Average Magnetization m (dimensionless)');  
title('Average Magnetization vs Temperature'); 
grid on; 

% Plot 2: Average Energy vs Temperature
figure(2);
plot(T_mc,e_bar,'b.', 'MarkerSize', 10);
xlabel('Temperature T (J)');  
ylabel('Average Energy E (J)');  
title('Average Energy vs Temperature'); 
grid on;

% Plot 3: Heat Capacity vs Temperature
figure(3);
plot(T_mc,cv_bar,'r.', 'MarkerSize', 10);
xlabel('Temperature T (J)'); 
ylabel('Heat capacity C (dimensionless)');  
title('Heat capacity vs Temperature'); 
grid on;

