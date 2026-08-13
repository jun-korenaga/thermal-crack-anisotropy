function d = load_numeric_table(filename,ncol)
% LOAD_NUMERIC_TABLE Load and validate an ASCII numeric table.

if ~exist(filename,'file')
  error('missing input file: %s',filename);
end
d = load(filename);
if ~isnumeric(d) || size(d,2)~=ncol || any(~isfinite(d(:)))
  error('invalid numeric table: %s',filename);
end
end
