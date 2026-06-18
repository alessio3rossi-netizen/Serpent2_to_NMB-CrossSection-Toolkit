clearvars

%% === Loop over all cases ===
rootDataDir = fullfile(pwd, "data");

fprintf("Current folder: %s\n", pwd);
fprintf("Looking for data in: %s\n", rootDataDir);

if ~exist(rootDataDir, "dir")
    error("Data folder does not exist: %s", rootDataDir);
end

libDirs = dir(rootDataDir);
libDirs = libDirs([libDirs.isdir] & ~startsWith({libDirs.name}, "."));

if isempty(libDirs)
    error("No library folders found inside: %s", rootDataDir);
end

for iLib = 1:numel(libDirs)
    library = string(libDirs(iLib).name);

    reacDirs = dir(fullfile(rootDataDir, library));
    reacDirs = reacDirs([reacDirs.isdir] & ~startsWith({reacDirs.name}, "."));

    for iReac = 1:numel(reacDirs)
        reactor = string(reacDirs(iReac).name);

        enrDirs = dir(fullfile(rootDataDir, library, reactor));
        enrDirs = enrDirs([enrDirs.isdir] & ~startsWith({enrDirs.name}, "."));

        for iEnr = 1:numel(enrDirs)
            enrichment = string(enrDirs(iEnr).name);

            %% === Case definition ===
            nameTag = library + "_" + reactor + "_" + enrichment;

            dataDir   = fullfile(rootDataDir, library, reactor, enrichment);
            resultDir = fullfile(pwd, "results", library, reactor, enrichment);
            if ~exist(resultDir, "dir"); mkdir(resultDir); end

            fprintf("\n=== Processing %s ===\n", nameTag);

            clear ZAI_BU Name_BU BU timeYears MDENS_BU ZAI_NMB NMB_names MDENS_NMB
            clear ZAI MAT_fuel_BURNUP MAT_fuel_MDENS MAT_UO2_BURNUP MAT_UO2_MDENS DAYS

            depFiles = dir(fullfile(dataDir, "*_dep.m"));

            if numel(depFiles) ~= 1
                error("Expected exactly one *_dep.m file in %s, found %d.", ...
                    dataDir, numel(depFiles));
            end

            depFile = fullfile(depFiles(1).folder, depFiles(1).name);

            outXlsx = fullfile(resultDir, nameTag + "_SerpentInventory_MDENS.xlsx");

            %% ===== Load required support modules =====
            run(fullfile("support_modules", "NMB_Isotopes.m"));  % must define NMB_isProduced

            if exist("ZAI2Name", "file") ~= 2
                error("ZAI2Name function was not found on the MATLAB path.");
            end

            %% ===== Run depletion file =====
            if ~exist(depFile, "file")
                error("Could not find depletion file: %s", depFile);
            end

            run(depFile);

            %% ===== Required common variables =====
            if ~exist("ZAI", "var")
                error("%s did not define ZAI.", depFiles(1).name);
            end

            if ~exist("DAYS", "var")
                error("%s did not define DAYS.", depFiles(1).name);
            end

            %% ===== Select material basis: fuel first, then UO2 =====
            if exist("MAT_fuel_BURNUP", "var") && exist("MAT_fuel_MDENS", "var")
                materialBasis = "fuel";
                BU_raw       = MAT_fuel_BURNUP;
                MDENS_full   = MAT_fuel_MDENS;

            elseif exist("MAT_UO2_BURNUP", "var") && exist("MAT_UO2_MDENS", "var")
                materialBasis = "UO2";
                BU_raw       = MAT_UO2_BURNUP;
                MDENS_full   = MAT_UO2_MDENS;

            else
                error(["%s did not define either fuel or UO2 depletion quantities.\n" + ...
                       "Expected either MAT_fuel_BURNUP and MAT_fuel_MDENS,\n" + ...
                       "or MAT_UO2_BURNUP and MAT_UO2_MDENS."], depFiles(1).name);
            end

            fprintf("Using material basis: MAT_%s_*\n", materialBasis);

            %% ===== Drop last two non-isotope entries =====
            ZAI_full = ZAI(:);
            nFull = numel(ZAI_full);

            if nFull < 3
                error("ZAI has only %d entries; cannot drop last two.", nFull);
            end

            isoIdx = 1:(nFull-2);

            ZAI_BU = ZAI_full(isoIdx);

            if size(MDENS_full, 1) < numel(isoIdx)
                error("MAT_%s_MDENS has %d rows but expected at least %d to match ZAI.", ...
                    materialBasis, size(MDENS_full,1), numel(isoIdx));
            end

            MDENS_BU = MDENS_full(isoIdx, :);

            BU = BU_raw(:);
            nBU = numel(BU);

            timeYears = DAYS(:) / 365.25;

            if size(MDENS_BU, 2) ~= nBU
                warning("MDENS cols (%d) != nBU (%d). Check depletion file.", ...
                    size(MDENS_BU,2), nBU);
            end

            if numel(timeYears) ~= nBU
                warning("numel(DAYS) (%d) != nBU (%d). Time headers may not align.", ...
                    numel(timeYears), nBU);
            end

            %% ===== Reorder MDENS into NMB isotope order =====
            Name_BU   = string(ZAI2Name(ZAI_BU));
            NMB_names = string(NMB_isProduced(:));

            [tf, locBU] = ismember(NMB_names, Name_BU);

            nNMB = numel(NMB_names);
            MDENS_NMB = nan(nNMB, size(MDENS_BU,2));
            ZAI_NMB   = nan(nNMB,1);

            MDENS_NMB(tf,:) = MDENS_BU(locBU(tf), :);
            ZAI_NMB(tf)     = ZAI_BU(locBU(tf));

            nMissing = sum(~tf);
            if nMissing > 0
                fprintf("Note: %d NMB isotopes were not found and will be NaN.\n", nMissing);
            end

            %% ===== Write Excel =====
            sheet = "inventory";

            writecell({char(nameTag) + " Serpent inventory: MAT_" + char(materialBasis) + "_MDENS (NMB order)"}, ...
                outXlsx, "Sheet", sheet, "Range", "A1");

            writecell({""}, outXlsx, "Sheet", sheet, "Range", "A2");
            writecell(cellstr("BU=" + string(BU.')), outXlsx, "Sheet", sheet, "Range", "B2");

            writecell({""}, outXlsx, "Sheet", sheet, "Range", "A3");
            writecell(cellstr("t[y]=" + string(timeYears.')), outXlsx, "Sheet", sheet, "Range", "B3");

            writecell(cellstr(NMB_names), outXlsx, "Sheet", sheet, "Range", "A4");
            writematrix(MDENS_NMB, outXlsx, "Sheet", sheet, "Range", "B4");

            %% ===== Save MAT =====
            save(fullfile(resultDir, nameTag + "_SerpentInventory.mat"), ...
                "materialBasis", ...
                "ZAI_BU", "Name_BU", "BU", "timeYears", "MDENS_BU", ...
                "ZAI_NMB", "NMB_names", "MDENS_NMB", "-v7.3");

            fprintf("\nInventory export complete:\n  Excel: %s\n  MAT  : %s\n", ...
                outXlsx, fullfile(resultDir, nameTag + "_SerpentInventory.mat"));
        end
    end
end