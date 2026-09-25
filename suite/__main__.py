# python -m suite
import sys
from . import core


def main():
    args = sys.argv[1:]
    if args and args[0] == "find" and len(args) > 1:
        for f in core.find(target=args[1]):
            print(f); print()
        return
    if args and args[0] == "dim" and len(args) > 1:
        for f in core.find(dimension=args[1]):
            print(f); print()
        return
    if args and args[0] == "report":
        from . import report
        print(report.run_report(target=args[1] if len(args) > 1 else None))
        return
    print("# q2sl suite (simulation)")
    print(f"# dimensions: {', '.join(core.dimensions())}")
    print(f"# targets   : {', '.join(core.targets())}")
    print("# modules:")
    for n, meta in core.catalog().items():
        print(f"#   {n:>22}  [{meta['dimension']}/{meta['target']}]  {meta['summary']}")
    print("\n# running every module once (defaults):\n")
    for f in core.find():
        print(f); print()
    print("# target one thing:  python -m suite find gate-identity")
    print("# or one dimension:  python -m suite dim syndrome")


if __name__ == "__main__":
    main()
