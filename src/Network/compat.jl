import Compat
using Compat: @warn, argmax, argmin, copyto!, findall, undef
using Markdown: @doc_str

if isdefined(Base, :NamedTuple)
    kwargs_pairs(kwargs::NamedTuple) = pairs(kwargs)
end
kwargs_pairs(kwargs) = kwargs
