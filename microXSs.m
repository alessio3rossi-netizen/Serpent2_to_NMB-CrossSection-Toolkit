clearvars

%% === Select cases to run ===
rootDataDir = fullfile(pwd, "data");

libraries = dir(rootDataDir);
libraries = libraries([libraries.isdir] & ~startsWith({libraries.name}, "."));

%% === For printing a common xlsx file ===
combinedXlsx = fullfile(pwd, "results", "ALL_microXS.xlsx");
combinedSheet = "microXS";
combinedRow = 1;

for iLib = 1:numel(libraries)
    library = string(libraries(iLib).name);

    reactorDirs = dir(fullfile(rootDataDir, library));
    reactorDirs = reactorDirs([reactorDirs.isdir] & ~startsWith({reactorDirs.name}, "."));

    for iReac = 1:numel(reactorDirs)
        reactor = string(reactorDirs(iReac).name);

        enrichDirs = dir(fullfile(rootDataDir, library, reactor));
        enrichDirs = enrichDirs([enrichDirs.isdir] & ~startsWith({enrichDirs.name}, "."));

        for iEnr = 1:numel(enrichDirs)
            enrichment = string(enrichDirs(iEnr).name);

            caseTag = library + "_" + reactor + "_" + enrichment;

            dataDir = fullfile(rootDataDir, library, reactor, enrichment);
            resultDir = fullfile(pwd, "results", library, reactor, enrichment);
            if ~exist(resultDir, "dir"); mkdir(resultDir); end

            fprintf("\n\n=== Running case: %s ===\n", caseTag);

            clear xsByBU xsActByBU allZAI allName masterSet avgT actXS validIdx BU BUv

            fprintf("Data directory: %s\n", dataDir);

            %% ---- your existing code starts here ----

matFile = fullfile("support_modules", library + ".mat");
if ~isfile(matFile)
    error("Decay matrix file not found: %s", matFile);
end
load(matFile); % Binary files containing the decay matrix of various libraries

run("support_modules\NMB_Isotopes.m"); % List of tabulated NMB isotopes
run("support_modules\nu_ORIGEN.m"); % neutron yield from ORIGEN
nu = nu_thermal; % !!! User input !!!

%% === Load burnup grid from _dep.m ===
depFiles = dir(fullfile(dataDir, "*_dep.m"));

if numel(depFiles) ~= 1
    error("Expected exactly one *_serp_dep.m file in %s, found %d.", ...
        dataDir, numel(depFiles));
end

depFile = fullfile(depFiles(1).folder, depFiles(1).name);
run(depFile);

if ~exist("MAT_fuel_BURNUP", "var")
    error("*_dep.m did not define MAT_fuel_BURNUP");
end
BU  = MAT_fuel_BURNUP(:); % BU values 
nBU = numel(BU); % Number of BU steps

%% === Find depmtx files and sort by bu_step ===
depmtxPattern = "*_depmtx_*_*_0.m"; % !!! User input !!!
L = dir(fullfile(dataDir, depmtxPattern));
if isempty(L)
    error("No files matched %s in %s", depmtxPattern, dataDir);
end

steps = nan(numel(L),1); % Preallocate an array for BU steps indices
for k = 1:numel(L) 
    tok = regexp(L(k).name, "_(\d+)_0\.m$", "tokens", "once"); 
    if ~isempty(tok)
        steps(k) = str2double(tok{1}); 
    end
end

keep = ~isnan(steps);
L = L(keep);
steps = steps(keep);

[steps, ord] = sort(steps); % sort if files are not in order
L = L(ord); 

fprintf("Found %d depmtx files (steps %d..%d).\n", numel(L), steps(1), steps(end));

%% === Storage ===
xsByBU     = cell(nBU, 1);   % FULL NMB table per burnup index
xsActByBU  = cell(nBU, 1);   % ACTINIDES ONLY table per burnup index

% Order-preserving master isotope lists
allZAI    = [];        % will be initialized from first valid table in NMB order
allName   = strings(0,1);
masterSet = containers.Map('KeyType','double','ValueType','logical'); % presence set

%% === Loop over depmtx files ===
for k = 1:numel(L)
    step  = steps(k);
    buIdx = step + 1;  % assumes step 0 corresponds to MAT_fuel_BURNUP(1)

    if buIdx < 1 || buIdx > nBU
        warning("Step %d maps to BU index %d (out of range). Skipping.", step, buIdx);
        continue;
    end

    filePath = fullfile(L(k).folder, L(k).name);
    fprintf("Processing step=%d (BU=%g): %s\n", step, BU(buIdx), L(k).name);

    run(filePath);

    if ~exist("A", "var")
        warning("File %s did not define variable A. Skipping.\n", L(k).name);
        continue;
    end

    result = extract_microXS(A, flx, ZAI, ZAI_lib, decayMatrix_lib, NMB_isProduced, nu);

    % --- Convert result.dataNMB -> table ---
    fields = cellstr(result.fields);
    D = result.dataNMB;

    ZAIcol  = double(D(:,1));
    Namecol = string(D(:,2));

    T = table(ZAIcol, Namecol, 'VariableNames', fields(1:2));
    for c = 3:numel(fields)
        T.(fields{c}) = double(D(:,c));
    end

    % Store FULL table
    xsByBU{buIdx} = T;

    % === Build/extend master isotope order WITHOUT sorting ===
    if isempty(allZAI)
        % Initialize from first valid table: keeps NMB order exactly
        allZAI  = T.ZAI(:);
        allName = T.Name(:);
        for ii = 1:numel(allZAI)
            masterSet(allZAI(ii)) = true;
        end
    else
        % Append only new isotopes (if any), keeping appearance order in T
        for ii = 1:height(T)
            zai = T.ZAI(ii);
            if ~isKey(masterSet, zai)
                allZAI(end+1,1)  = zai;         %#ok<SAGROW>
                allName(end+1,1) = T.Name(ii);  %#ok<SAGROW>
                masterSet(zai) = true;
            end
        end
    end

    % --- Actinides-only subset ---
    Znum = floor(T.ZAI / 10000);
    isAct = (Znum >= 89) & (Znum <= 103);
    xsActByBU{buIdx} = T(isAct, :);

    clear A
end

%% === Burnup-weighted averages (ALL isotopes, in ORIGINAL order) ===
validIdx = find(~cellfun(@isempty, xsByBU));
if numel(validIdx) < 2
    error("Need at least two burnup points with data to compute weighted averages.");
end

BUv = BU(validIdx);

% Step-average weighting: XS_k applies over interval [BU_k, BU_{k+1}]
dBU = diff(BUv);
if any(dBU <= 0)
    error("MAT_fuel_BURNUP must be strictly increasing over used points.");
end

xsCols = fields(3:end);

% Output table in the preserved order
avgT = table(allZAI(:), allName(:), 'VariableNames', {'ZAI','Name'});

for c = 1:numel(xsCols)
    col = xsCols{c};
    avgVals = nan(height(avgT), 1);

    for i = 1:height(avgT)
        zai = avgT.ZAI(i);

        xs = nan(numel(validIdx), 1);
        for j = 1:numel(validIdx)
            idx = validIdx(j);
            T = xsByBU{idx};
            hit = find(T.ZAI == zai, 1);
            if ~isempty(hit)
                xs(j) = T.(col)(hit);
            end
        end

        xs_int = xs(1:end-1);
        m = isfinite(xs_int);

        if any(m)
            avgVals(i) = sum(xs_int(m) .* dBU(m)) / sum(dBU(m));
        else
            avgVals(i) = NaN;
        end
    end

    avgT.(col) = avgVals;
end

%% === Save outputs ===
save(fullfile(resultDir, caseTag + "_microXS_byBU_FULL_and_ACTINIDES.mat"), ...
    "xsByBU", "xsActByBU", "BU", "validIdx", "-v7.3");

save(fullfile(resultDir, caseTag + "_microXS_weightedAvg_ALL_isotopes_ORDERED.mat"), ...
    "avgT", "BU", "validIdx");

fprintf("\nSaved (order-preserving):\n  - Burnup-dependent FULL + actinides-only\n  - Burnup-weighted average ALL isotopes (same order as first table)\nResults in: %s\n", resultDir);


%% ===== Export microXS to Excel =====
% Required in workspace: avgT, xsActByBU, BU, validIdx
%
% avgT columns: ZAI, Name, and reaction columns: NF,F,A,N2N,AEx,N2NEx,N3N,NAlpha,NP
% xsActByBU{idx}: table for that BU index (actinides only) with same columns
% BU: burnup vector, validIdx: indices actually processed

% ---- Define the <name> tag (used in headers) ----
nameTag = caseTag;      % <-- set your <name> here (example from your screenshots)

% ---- Output path ----
if ~exist(resultDir, "dir"); mkdir(resultDir); end
outXlsx = fullfile(resultDir, nameTag + "_microXS.xlsx");

% ---- Reaction columns (must match avgT variable names) ----
rxns = ["NF","F","A","N2N","AEx","N2NEx","N3N","NAlpha","NP"];

% ---- Burnup points used ----
validIdx = validIdx(:);
BUv = BU(validIdx);
nBUv = numel(BUv);

% ---- Actinide isotope order (reference) ----
firstAct = find(~cellfun(@isempty, xsActByBU), 1, "first");
if isempty(firstAct)
    error("xsActByBU is empty everywhere. Cannot export actinide blocks.");
end
Tact0  = xsActByBU{firstAct};
actZAI  = Tact0.ZAI(:);
actName = string(Tact0.Name(:));
nAct = numel(actZAI);

% ---- Build burnup-dependent matrices for each reaction (actinides only) ----
actXS = struct();
for r = 1:numel(rxns)
    actXS.(rxns(r)) = nan(nAct, nBUv);
end

for j = 1:nBUv
    idx = validIdx(j);
    T = xsActByBU{idx};
    if isempty(T); continue; end

    [tf, loc] = ismember(actZAI, T.ZAI);
    for r = 1:numel(rxns)
        col = rxns(r);
        vals = nan(nAct,1);
        vals(tf) = T.(col)(loc(tf));
        actXS.(col)(:, j) = vals;
    end
end

%% ===== Write to Excel =====
sheet = "microXS";

% (A) AVERAGE TABLE (all isotopes, SAME ORDER as avgT)
% Layout:
%   A1: <name>
%   A2: blank, B2..: X_<name>_<RXN>
%   A3..: isotope names, B3..: values

% Title
writecell({char("X_" + nameTag)}, outXlsx, "Sheet", sheet, "Range", "A1");

% Headers (no "Isotope" tag)
avgHeaders = "X_" + nameTag + "_" + rxns;  % e.g. X_PWR34_40_NF
writecell({""}, outXlsx, "Sheet", sheet, "Range", "A2");
writecell(cellstr(avgHeaders), outXlsx, "Sheet", sheet, "Range", "B2");

% Isotope names column, with artificial FP at the end
isoLabelAll = [string(avgT.Name); "FP"];

writecell(cellstr(isoLabelAll), outXlsx, "Sheet", sheet, "Range", "A3");

% Numeric matrix, with artificial FP row of zeros at the end, for format
avgNum = nan(height(avgT), numel(rxns));
for r = 1:numel(rxns)
    avgNum(:, r) = avgT.(rxns(r));
end

avgNum = [avgNum; zeros(1, numel(rxns))];

writematrix(avgNum, outXlsx, "Sheet", sheet, "Range", "B3");

% (B) BU-DEPENDENT BLOCKS (actinides only), side-by-side
% Each block layout:
%   <col><row>     : X_<name>_<RXN>_BD
%   <col><row+1>   : BD | BU values ...
%   <col><row+2..> : actinide names | XS(BU) ...

startRow = 5 + height(avgT) + 1; % leave a little space under averages, one extra line for FPs (will be empty)

blockWidth = 1 + nBUv;  % first col = names, next nBUv cols = BU points
gap = 2;                % blank columns between blocks

for r = 1:numel(rxns)
    rxn = rxns(r);

    startCol = 1 + (r-1) * (blockWidth + gap);

    % Title: X_<name>_<RXN>_BD
    writecell({char("X_" + nameTag + "_" + rxn + "_BD")}, outXlsx, ...
        "Sheet", sheet, "Range", excelColumn(startCol) + string(startRow));

    % Header row: BD + burnup values
    hdrRow = startRow + 1;
    writecell({char("BD")}, outXlsx, "Sheet", sheet, ...
        "Range", excelColumn(startCol) + string(hdrRow));
    writematrix(BUv.', outXlsx, "Sheet", sheet, ...
        "Range", excelColumn(startCol+1) + string(hdrRow));

    % Data rows: actinide names + XS matrix
    dataRow = startRow + 2;
    writecell(cellstr(actName), outXlsx, "Sheet", sheet, ...
        "Range", excelColumn(startCol) + string(dataRow));

    writematrix(actXS.(rxn), outXlsx, "Sheet", sheet, ...
        "Range", excelColumn(startCol+1) + string(dataRow));
end

%% ===== Common xlxs file =====

% --- Averaged XS block ---
writecell({char(nameTag)}, combinedXlsx, ...
    "Sheet", combinedSheet, ...
    "Range", "A" + string(combinedRow));

writecell({""}, combinedXlsx, ...
    "Sheet", combinedSheet, ...
    "Range", "A" + string(combinedRow + 1));

writecell(cellstr(avgHeaders), combinedXlsx, ...
    "Sheet", combinedSheet, ...
    "Range", "B" + string(combinedRow + 1));

writecell(cellstr(isoLabelAll), combinedXlsx, ...
    "Sheet", combinedSheet, ...
    "Range", "A" + string(combinedRow + 2));

writematrix(avgNum, combinedXlsx, ...
    "Sheet", combinedSheet, ...
    "Range", "B" + string(combinedRow + 2));

% Move below averaged block
combinedRow = combinedRow + numel(isoLabelAll) + 4;


% --- BU-dependent actinide blocks, side-by-side ---
bdStartRow = combinedRow;

for r = 1:numel(rxns)
    rxn = rxns(r);

    startCol = 1 + (r-1) * (blockWidth + gap);

    % Title
    writecell({char("X_" + nameTag + "_" + rxn + "_BD")}, combinedXlsx, ...
        "Sheet", combinedSheet, ...
        "Range", excelColumn(startCol) + string(bdStartRow));

    % Header row: BD + BU values
    hdrRow = bdStartRow + 1;

    writecell({char("BD")}, combinedXlsx, ...
        "Sheet", combinedSheet, ...
        "Range", excelColumn(startCol) + string(hdrRow));

    writematrix(BUv.', combinedXlsx, ...
        "Sheet", combinedSheet, ...
        "Range", excelColumn(startCol + 1) + string(hdrRow));

    % Data rows
    dataRow = bdStartRow + 2;

    writecell(cellstr(actName), combinedXlsx, ...
        "Sheet", combinedSheet, ...
        "Range", excelColumn(startCol) + string(dataRow));

    writematrix(actXS.(rxn), combinedXlsx, ...
        "Sheet", combinedSheet, ...
        "Range", excelColumn(startCol + 1) + string(dataRow));
end

% Move below entire side-by-side BD section
combinedRow = bdStartRow + nAct + 5;
%% ---- existing code ends here ----
        end
    end
end

%% ===== Helper: convert column index -> Excel letters =====
function letters = excelColumn(colIdx)
    letters = "";
    while colIdx > 0
        rem = mod(colIdx - 1, 26);
        letters = char(65 + rem) + letters;
        colIdx = floor((colIdx - 1) / 26);
    end
end

