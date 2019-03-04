if __name__ == "__main__":

    from dials_scratch.jmp.profile_viewer.extract_strong import extract_strong
    from dials_scratch.jmp.profile_viewer.extract_strong import generate_profile
    from dials_scratch.jmp.profile_viewer.extract_strong import show_profile
    from dials_scratch.jmp.profile_viewer.extract_strong import crystal_parameters
    import argparse

    parser = argparse.ArgumentParser(description="Extract and view profiles")

    parser.add_argument(
        "--directory",
        dest="directory",
        default=None,
        help="The directory which must contain integrated_experiments.json and shoeboxes*.pickle files",
    )

    parser.add_argument(
        "--start",
        dest="start",
        choices=["extract", "generate", "display", "decompose"],
        default="extract",
        help="If data has been generated then start later",
    )

    parser.add_argument(
        "--sample",
        dest="sample",
        type=int,
        default=None,
        help="Choose a random sample of reflections",
    )

    parser.add_argument(
        "--grid_size",
        dest="grid_size",
        type=int,
        default=15,
        help="The grid size is 2N+1",
    )

    parser.add_argument(
        "--grid_range",
        dest="grid_range",
        type=float,
        default=0.25,
        help="The grid range (+-HKL)",
    )

    args = parser.parse_args()

    if args.start == "extract":
        extract_strong(args.directory)
    if args.start in ["extract", "generate"]:
        generate_profile(args.sample, args.grid_size, args.grid_range)
    if args.start in ["extract", "generate", "display"]:
        show_profile(args.grid_range)
    crystal_parameters()
