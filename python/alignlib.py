# ignore warnings for double-defined templates
# like "to-Python converter for std::vector<double, std::allocator<double> > already registered; second conversion method ignored"
import warnings 
warnings.filterwarnings("ignore")
from _alignlib import *
