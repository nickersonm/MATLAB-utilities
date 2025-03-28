function package = packageDependencies(scriptName)

file = dir(scriptName);

% Get file dependencies
deps = matlab.codetools.requiredFilesAndProducts(file.folder + "/" + file.name);
deps = unique(deps(:));

% Assemble package
[~, package] = fileparts(scriptName);
<<<<<<< HEAD
package = "'" + string(package) + "'_and_dependencies.zip";
=======
package = string(package) + "_and_dependencies.zip";
>>>>>>> 54d215d43c3f26fd59fd89a403fd874e2bb21b97
zip(package, deps);

fprintf("\nPackaged %i dependencies of '%s' into '%s'\n", numel(deps), scriptName, package);

end
