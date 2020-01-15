project: ALaDyn
src_dir: ./src
output_dir: ./html/doc_ford
exclude_dir: ./src/depot
include: /usr/include
         /usr/local/include
         /usr/include/mpi
summary: A High-Accuracy PIC Code for the Maxwell-Vlasov Equations
author: The ALaDyn Collaboration
author_description:
project_github: https://github.com/ALaDyn/ALaDyn
project_download: https://github.com/ALaDyn/ALaDyn/releases/latest
github: https://github.com/ALaDyn
email: 
media: ./media
docmark: !
predocmark: >
docmark_alt: #
predocmark_alt: <
display: public
	protected
	private
source: true
graph: true
search: true
warn: false
license: gfdl
version: {!./docs/VERSION.md!}
page_dir: ./docs/pages
favicon: ./media/site_icon.png
coloured_edges: true
print_creation_date: true
creation_date: %Y-%m-%d

{!./docs/INTRO.md!}

[**General code guide**](|url|/page/index.html)  

[Papers published by the ALaDyn Collaboration](https://aladyn.github.io/Papers/)

[Code description](|url|/page/DESCRIPTION.html)  

[Input guide](|url|/page/NAMELIST_GUIDE.html)

@Bug
ALaDyn code is fully maintained and tested.
Anyway, if while using it you find any bug affecting the code or if you simply
have questions, please feel free to open an issue in the [GitHub page](https://github.com/ALaDyn/ALaDyn/issues).
@endbug

Copyright on the code is by the ALaDyn Collaboration.
