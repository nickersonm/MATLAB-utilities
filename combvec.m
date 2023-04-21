% Find all possible combinations of input vectors
% Michael Nickerson 2022-03-22
function varargin = combvec(varargin)
  [varargin{:}] = ndgrid(varargin{:});
  varargin = reshape(cell2mat(varargin),[],numel(varargin));
end
