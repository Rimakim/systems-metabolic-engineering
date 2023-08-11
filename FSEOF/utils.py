from argparse import ArgumentParser

def argument_parser(version=None):
    def boolean_string(s):
        if s not in {'False', 'True'}:
            raise ValueError('Not a valid boolean string')
        return s == 'True'
    
    parser = ArgumentParser()
    parser.add_argument('-i', '--model_dir', 
                        required=True, help='Input model directory')
    parser.add_argument('-o', '--output_dir', 
                        required=True, help='Output directory')
    parser.add_argument('-b', '--biomass_reaction', required=False,
                        default='None', help='BIOMASS reaction ID')
    parser.add_argument('-t', '--objective_reaction', required=True, 
                        help='Target chemical production reaction ID')
    parser.add_argument('-max_step', '--max_step', required=False, type=int,
                        default=10, 
                        help='Maximum number of knockout targets')
    parser.add_argument('-cpu', '--cpu_number', required=False, type=int,
                        default=8, 
                        help='Number of cpu thread to use')
    parser.add_argument('-oxygen', '--oxygen_condition', required=False, 
                        default='aerobic', 
                        choices=['aerobic', 'anaerobic'],
                        help='Air condition')
    return parser
