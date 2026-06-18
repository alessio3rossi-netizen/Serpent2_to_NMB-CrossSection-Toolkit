function result = extract_microXS(A, flx, ZAI, ZAI_lib, decayMatrix_lib, NMB_isProduced, nu)
%--------------------------------------------------------------------------
% extract_microXS 
%--------------------------------------------------------------------------

% ====== (1) Decay struct input ======
idx1 = ismember(ZAI, ZAI_lib); 
idx2 = ismember(ZAI_lib, ZAI);

NameSerpent = string(ZAI2Name(ZAI));

decayMatrix = zeros(length(ZAI));
decayMatrix(idx1, idx1) = decayMatrix_lib(idx2, idx2);

% ====== (2) Serpent matrix and flux ======
burnupMatrix = A;

irrMatrix = burnupMatrix - decayMatrix;
save("debug", "irrMatrix", "burnupMatrix", "decayMatrix", "ZAI", "ZAI_lib");

fuelFlux = flx;

% ====== (3) Correct NMB ordering ======
% tf = match exists; loc = index in NameIsProduced corresponding to each NMB entry
[tf, loc] = ismember(NMB_isProduced, NameSerpent);

% ZAI and Names in correct NMB order
ZAI_NMB = ZAI(loc);
Name_NMB = NameSerpent(loc);

% ====== (4) Initialize reaction XS arrays ======
nIso = length(ZAI);
fiss = zeros(nIso,1);
ngamma = zeros(nIso,1);
ngammaEX = zeros(nIso,1);
np = zeros(nIso,1);
nalpha = zeros(nIso,1);
n2n = zeros(nIso,1);
n2nEX = zeros(nIso,1);
n3n = zeros(nIso,1);

% ====== (5) Compute microscopic cross sections ======

for j = 1:nIso
    for i = 1:nIso
        if irrMatrix(i,j) ~= 0
            ZAIdiff = ZAI(j) - ZAI(i);

            % Fission (empirical)
            if sum(irrMatrix(:,j) ~= 0) > 100 && j - i < 50 

                fiss(j) = fiss(j) - irrMatrix(i,j);
            end

            % Reaction pattern recognition
            if ZAIdiff == -10
                ngamma(j) = ngamma(j) + irrMatrix(i,j);     % (n,gamma)
            elseif ismember(ZAIdiff,[-12, -11])
                ngammaEX(j) = ngammaEX(j) + irrMatrix(i,j);      % (n,gamma)EX
            elseif ZAIdiff == 10000
                np(j) = np(j) + irrMatrix(i,j);             % (n,p)
            elseif ZAIdiff == 20030
                nalpha(j) = nalpha(j) + irrMatrix(i,j);     % (n,α)
            elseif ZAIdiff == 10
                n2n(j) = n2n(j) + irrMatrix(i,j);           % (n,2n)
            elseif ismember(ZAIdiff,[8, 9])
                n2nEX(j) = n2nEX(j) + irrMatrix(i,j);           % (n,2n)EX
            elseif ZAIdiff == 20
                n3n(j) = n3n(j) + irrMatrix(i,j);           % (n,3n)
            end
        end
    end

    % Convert to barns
    fiss(j)   = (fiss(j)   / fuelFlux) * 1e24;
    ngamma(j) = (ngamma(j) / fuelFlux) * 1e24;
    ngammaEX(j) = (ngammaEX(j) / fuelFlux) * 1e24;
    np(j)     = (np(j)     / fuelFlux) * 1e24;
    nalpha(j) = (nalpha(j) / fuelFlux) * 1e24;
    n2n(j)    = (n2n(j)    / fuelFlux) * 1e24;
    n2nEX(j)    = (n2nEX(j)    / fuelFlux) * 1e24;
    n3n(j)    = (n3n(j)    / fuelFlux) * 1e24;
end

% ====== (6) Map XS arrays into NMB order ======
%totNMB    = tot(loc);
fissNMB   = fiss(loc);
nfissNMB = fissNMB;

nfissNMB(1:length(nu)) = fissNMB(1:length(nu)).*nu;

ngammaNMB = ngamma(loc);
ngammaEXNMB = ngammaEX(loc);
npNMB     = np(loc);
nalphaNMB = nalpha(loc);
n2nNMB    = n2n(loc);
n2nEXNMB    = n2nEX(loc);
n3nNMB    = n3n(loc);

% ====== (7) Collect outputs ======
fields = ["ZAI", "Name","NF","F", "A", "N2N", "AEx", "N2NEx", "N3N", "NAlpha", "NP"];

dataNMB = [ZAI_NMB, Name_NMB, nfissNMB, fissNMB, ngammaNMB, n2nNMB, ngammaEXNMB, n2nEXNMB, n3nNMB, nalphaNMB, npNMB];

% Store results
result.dataNMB   = dataNMB;
result.fields    = fields;
result.tableNMB = table( ...
    ZAI_NMB, Name_NMB, nfissNMB, fissNMB, ngammaNMB, n2nNMB, ngammaEXNMB, ...
    n2nEXNMB, n3nNMB, nalphaNMB, npNMB, ...
    'VariableNames', cellstr(fields));
%result.decayMatrix = decayMatrix;
%result.irrMatrix = irrMatrix;
%result.burnMatrix = burnupMatrix;

end