from argparse import ArgumentParser

def argument_parser(version=None):
    def boolean_string(s):
        if s not in {'False', 'True'}:
            raise ValueError('Not a valid boolean string')
        return s == 'True'
    
    parser = ArgumentParser()
    parser.add_argument('-m', '--model_dir', 
                        required=True, help='Input model directory')
    parser.add_argument('-b', '--biomass_lb', required=False, type=float,
                        default = 0.1,
                        help='minimal biomass production')
    parser.add_argument('-o', '--output_dir', required=True, 
                        help='output directory name')
    parser.add_argument('-t', '--objective_reaction', required=True, 
                        help='Target chemical production reaction ID')
    parser.add_argument('-max_ko', '--max_knockout', required=False, type=int,
                        default=1, 
                        help='Maximum number of knockout targets')
    parser.add_argument('-cpu', '--cpu_num', required=False, type=int,
                        default=8, 
                        help='Number of cpu thread to use')
    parser.add_argument('-oxygen', '--oxygen_condition', required=False, 
                        default='aerobic', 
                        choices=['aerobic', 'anaerobic'],
                        help='Air condition')
    parser.add_argument('-db', '--default_bound', required=False, type=float,
                        default=1000, 
                        help='bounds of dual variables')
    parser.add_argument('-it', '--iter_num', required='False', type=int,
                       default=1,
                       help = 'iteration number')
   
    return parser
