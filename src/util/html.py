#!/usr/bin/env python3

def print_html_header(fp_html):
    print("<!DOCTYPE html>", file=fp_html)
    print("<html>", file=fp_html)
    print("<head>", file=fp_html)
    print("<style>", file=fp_html)
    print("body {", file=fp_html)
    print("  font-family: \"Helvetica Neue\", \"Helvetica\", \"Arial\" ", file=fp_html)
    print("}", file=fp_html)
    print("</style>", file=fp_html)
    print("</head>", file=fp_html)
    print("<body>", file=fp_html)

def print_html_footer(fp_html):
    print("<br>", file=fp_html)    
    print("</body>", file=fp_html)
    print("</html>", file=fp_html)
        

