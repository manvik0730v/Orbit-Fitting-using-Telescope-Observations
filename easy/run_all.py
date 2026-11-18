"""Run the experiments, figures, slides, and checks from the repository root."""
from .experiments import run
from .plots import make_plots
from .slides import create
from .extras import run_extras
def main():
    """Generate all primary deliverables."""
    out,d,tau,u,fits,hard=run(); make_plots(); create(); run_extras(out,d,tau,u,fits)
if __name__=='__main__': main()
