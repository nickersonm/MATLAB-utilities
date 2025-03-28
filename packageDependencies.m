function package = packageDependencies(scriptName)

file = dir(scriptName);

% Get file dependencies
deps = matlab.codetools.requiredFilesAndProducts(file.folder + "/" + file.name);
deps = unique(deps(:));

% Assemble package
[~, package] = fileparts(scriptName);
package = "'" + string(package) + "'_and_dependencies.zip";
zip(package, deps);

fprintf("\nPackaged %i dependencies of '%s' into '%s'\n", numel(deps), scriptName, package);

end
