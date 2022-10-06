% Izhikevich model used up to Padova paper submission

function [spikeTrains,s_time,post,delays]=Izhikevich_network_v1(minutes,M,D,Ne,Ni)
% spnet.m: Spiking network with axonal conduction delays and STDP
% Created by Eugene M.Izhikevich.                February 3, 2004
% Modified to allow arbitrary delay distributions.  April 16,2008
N = Ne+Ni;

a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=10;                 % maximal synaptic strength

seconds = minutes*60;
totTime = 1000*seconds; % in ms

% post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
% Take special care not to have multiple connections between neurons
delays = cell(N,D);
for i=1:Ne
    p=randperm(N);
    post(i,:)= p(1:M);
    for j=1:M
        delays{i, ceil(D*rand)}(end+1) = j;  % Assign random exc delays
    end;
end;
for i=Ne+1:N
    p=randperm(Ne); % No inh to inh connections
    post(i,:)=p(1:M);
    delays{i,1}=1:M;                    % all inh delays are 1 ms.
end;

s=[6*ones(Ne,M);-5*ones(Ni,M)];         % synaptic weights
sd=zeros(N,M);                          % their derivatives
s_time=zeros(N,M,seconds); 

% Make links at postsynaptic targets to the presynaptic weights
pre = cell(N,1);
aux = cell(N,1);
for i=1:Ne
    for j=1:D
        for k=1:length(delays{i,j})
            pre{post(i, delays{i, j}(k))}(end+1) = N*(delays{i, j}(k)-1)+i;
            aux{post(i, delays{i, j}(k))}(end+1) = N*(D-1-j)+i; % takes into account delay
        end;
    end;
end;
  

STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings
spikeTrains = zeros(N, totTime);
numFired = zeros(1, totTime);

%for sec=1:60*60*24                      % simulation of 1 day
for minu = 1:minutes
    disp([num2str(minu), ' minutes'])
    for sec=1:60                            % simulation of 1 min
        sec_tot = 60*(minu-1)+sec;
      for t=1:1000                          % simulation of 1 sec
        t_tot_global = ((minu-1)*60 + sec-1)*1000 + t;
        I=zeros(N,1);        
        I(ceil(N*rand))=20;                 % random thalamic input 
        fired = find(v>=30);                % indices of fired neurons
        numFired(t_tot_global) = numel(fired);
        v(fired)=-65;  
        u(fired)=u(fired)+d(fired);
        STDP(fired,t+D)=0.1;
        for k=1:length(fired)
          sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
        end;
        firings=[firings;t*ones(length(fired),1),fired];
        k=size(firings,1);
        while firings(k,1)>t-D
          del=delays{firings(k,2),t-firings(k,1)+1};
          ind = post(firings(k,2),del);
          I(ind)=I(ind)+s(firings(k,2), del)';
          sd(firings(k,2),del)=sd(firings(k,2),del)-1.2*STDP(ind,t+D)';
          k=k-1;
        end;
        v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
        v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
        u=u+a.*(0.2*v-u);                   % step is 0.5 ms
        STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
      end;
      %plot(firings(:,1),firings(:,2),'.');
      %axis([0 1000 0 N]); drawnow;
      for tmpT = 1:1000
        t_tot = ((minu-1)*60+sec-1)*1000 + tmpT;
        fIdx = find(firings(:,1) == tmpT);
        neuId = firings(fIdx,2);
        spikeTrains(neuId,t_tot) = 1;
      end
      STDP(:,1:D+1)=STDP(:,1001:1001+D);
      ind = find(firings(:,1) > 1001-D);
      firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
      % synaptic weight updated every second
      s(1:Ne,:)=max(0,min(sm,0.01+s(1:Ne,:)+sd(1:Ne,:)));
      s_time(:,:,sec_tot)=s;
      sd=0.9*sd;
    end;
end

end