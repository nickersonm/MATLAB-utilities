% Get file dependencies
deps = [];

for file = dir("./*.m")'
    deps = [deps, matlab.codetools.requiredFilesAndProducts(file.folder + "/" + file.name)];
end

deps = unique(deps(:));

writecell(deps, "deps.txt");
disp(deps);
