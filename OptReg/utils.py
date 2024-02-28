from argparse import ArgumentParser

def argument_parser(version=None):
    def boolean_string(s):
        if s not in {'False', 'True'}:
            raise ValueError('Not a valid boolean string')
        return s == 'True'
    
    parser = ArgumentParser()
    parser.add_argument('-m', '--model_path', 
                        required=True, help='Input model path')
    
    parser.add_argument('-b', '--biomass_lb', required=False, type=float,
                        default = 0.1,
                        help='minimal biomass production')
    
    parser.add_argument('-o', '--output_path', required=True, 
                        help='Input output path')
    
    parser.add_argument('-t', '--target_rxn', required=True, 
                        help='Target chemical production reaction ID')
    
    parser.add_argument('-max_m', '--max_manipulation', required=False, type=int,
                        default=1, 
                        help='Maximum number of manipulation targets')
    
    parser.add_argument('-str', '--regulation_strength', required=False, type=int,
                        default=0.5, 
                        help='Flux regulation strength')
    
    parser.add_argument('-e', '--eps', required=False, type=int,
                        default=0.0005, 
                        help='Maximum number of knockout targets')
    
    parser.add_argument('-oxygen', '--oxygen_condition', required=False, 
                        default='aerobic', 
                        choices=['aerobic', 'anaerobic'],
                        help='Air condition')

    
    parser.add_argument('-cpu', '--cpu_num', required=False, type=int,
                        default=8, 
                        help='Number of cpu thread to use')
    
    parser.add_argument('-db', '--default_bound', required=False, type=float,
                        default=1000, 
                        help='bounds of dual variables')
    parser.add_argument('-it', '--iter_num', required=False, type=int,
                       default=1,
                       help = 'iteration number')
   
    return parser
