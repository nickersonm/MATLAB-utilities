% Simple fallback for using GPU arrays
function x = gpuArrayTry(x)
    % Convert to gpuArrayTry if available
    if gpuDeviceCount > 0
        x = gpuArray(x);
    end
end
