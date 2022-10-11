from dataclasses import dataclass

@dataclass(frozen=True)
class HelpItem:
    """Class for storing functions and associated help"""
    name: str
    function: str
    help_function: str
    help: str
    desc: str
    
def setup_subparsers(parser, subcommand_data):
    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True
    for subcommand_item in subcommand_data:
        subparser = subparsers.add_parser( subcommand_item.name,
                        description=subcommand_item.desc, help=subcommand_item.help)
        subparser = subcommand_item.help_function(subparser)
        subparser.set_defaults(func=subcommand_item.function)
    return parser
 
