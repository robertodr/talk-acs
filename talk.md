name: inverse
layout: true
class: center, middle, inverse
---
## Reusable Components for Quantum Chemistry Software

.author[Roberto Di Remigio]

.institution[UiT - The Arctic University of Norway]

.institution[Department of Chemistry - Virginia Tech]

.date[1 April 2019, Orlando]

.footnote[[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) licensed.
Browse slides at [http://tinyurl.com/talk-orlando](http://tinyurl.com/talk-orlando)]

???

I will present very recent efforts in designing and implementing a library for
input parsing that can be reused by different programs with minimal effort.
The library allows combination of different input formats and enforces sane
defaults, type checking, and documentation.
We achieve this with a strict separation between tokenization and validation.

---

layout: false

## Software Grows Complex Over Time.red[<sup>1</sup>]

_"The complexity of software is an essential property, not an accidental one."_.red[<sup>2</sup>]

<p style="text-align:center;"><img src="images/complex-machine.jpg" style="width: 30%"></p>
<p style="clear: both;">

.footnote-cite[.red[<sup>1</sup>] Hatton, L. _IEEE Trans. Software Eng._ (2014), **40**, 450

.red[<sup>2</sup>] Brooks, F. P., Jr. [_The Mythical Man-Month: Essays on Software Engineering_](http://worrydream.com/refs/Brooks-NoSilverBullet.pdf), Addison-Wesley Professional, 1995]

???

- We start with a simple idea and overtime software grows complex
- It is inherent characteristic of software to change and grow over time
- There is no silver bullet for managing complexity. However, there are some bad practices that we can avoid.

---

## Composition

- Build complex behavior from simple components
- We can reason about the components and the composite
- Composition is key to managing complexity
- Modularity does not imply simplicity, but is enabled by it

<img src="images/knit_vs_lego.jpg" style="width: 100%;"/>


.footnote-cite[Slide taken from [Complexity in software development by Jonas Juselius](https://github.com/scisoft/complexity)]

---
layout: false

## Extensible
### New Challenges! New Methods! New Software?!?

<p style="text-align:center;"><img src="images/yoda-future-once.gif" style="width: 100%;"/></p>
<p style="clear: both;">

???

- New modelling challenges require new (or improved) methods and new (or
  improved) software implementations.
- Software cannot be rewritten from scratch every time, we need to have
  _extensible_ architectures

---
layout: false

## Modular
### Do One Thing, Do It Well
<p style="text-align:center;"><img src="images/lego-once.gif" style="width: 100%"></p>
<p style="clear: both;">

???

- Fast-paced development cycles can only occur within a clear software structure
- Delimit capabilities, concerns, results
- Modularity can help create such a structure
- A LEGO-like structure can enhance our capabilities to reason about complex software workflows.
- Assemble, de-assemble, re-assemble is easier with standalone modules.

---
layout: false

## Sustainable
### Software Reuse for a Better Life

<p style="text-align:center;"><img src="images/reinvent-wheel.gif" style="width: 100%;"/></p>
<p style="clear: both;">

.footnote-cite[Source: [IBM Big Data & Analytics Hub](http://www.ibmbigdatahub.com/blog/dont-reinvent-wheel-increase-productivity-strategic-reuse)]

???

- Reinventing the wheel is bad
- Build small, well-documented software components
- Reuse components in different contexts

---

## The Problem(s) with Input Parsing

- Proliferation of formats
- Documentation quickly outdated or incomplete
- Hard to compose

---
layout: false
class: split-50-50

## parselglossy: generic input parsing, speaking in tongues

Separate _syntax_ and _semantics_

.column[
- [github.com/dev-cafe/parselglossy](https://github.com/dev-cafe/parselglossy)
- [parselglossy.readthedocs.io](https://parselglossy.readthedocs.io)
- `pip install parselglossy`
]
.column[
<p style="text-align:center;"><img src="images/parseltongue.jpg" style="width: 35%;"/></p>
<p style="clear: both;">
]

```
parselglossy --help

Usage: parselglossy [OPTIONS] COMMAND [ARGS]...

  Console script for parselglossy.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  doc       Generate documentation from validation template.
  lex       Run lexer and dump intermediate representation to JSON.
  parse     Parse input file.
  validate  Validate intermediate representation into final representation.
```

---
layout: false

## parselglossy

<p style="text-align:center;"><img src="images/structure.png" style="width: 100%;"/></p>
<p style="clear: both;">

---
layout: false

## Tokenization: how _any_ correct input _looks_ 

---
layout: false

## Validation: what correct input _must contain_

<p style="text-align:center;"><img src="images/template.png" style="width: 100%;"/></p>
<p style="clear: both;">

---
layout: false

<p style="text-align:center;"><img src="images/validation.png" style="width: 100%;"/></p>
<p style="clear: both;">

---
layout: false

## Bonus: autogenerate documentation

```
parselglossy doc template.yml --outfile input.rst
```

<p style="text-align:right;"><img src="images/documentation.png" style="width: 60%;"/></p>
<p style="clear: both;">

---
layout: false
class: split-40-60

## Acknowledgements

.column[
<p style="text-align:left;"><img src="images/radovan-bast.jpg" style="width: 80%"></p>
<p style="clear: both;">
]
.column[
__Radovan Bast__ .cite[UiT -- The Arctic University of Norway]
]

.column[
__Stig Rune Jensen__ .cite[Hylleraas Centre, UiT]
]
.column[
<p style="text-align:left;"><img src="images/stig-rune-jensen.jpeg" style="width: 10%"></p>
<p style="clear: both;">
]

.column[
__Jonas Juselius__ .cite[]
]
.column[
<p style="text-align:left;"><img src="images/jonas.jpg" style="width: 20%"></p>
<p style="clear: both;">
]

<p style="text-align:center;"><img src="images/NFR-logo.png" style="width: 80%"></p>
<p style="clear: both;">

---
layout: false

## Commercial

All examples online [github.com/dev-cafe/cmake-cookbook](https://github.com/dev-cafe/cmake-cookbook)

<p style="text-align:center;"><img src="images/cmake-cookbook-cover.png" style="width: 50%"></p>
<p style="clear: both;">

---
name: last-page
template: inverse

# Thanks for your attention!

Slideshow created using [remark] and served using [cicero]

Slides available on [GitHub](https://github.com/robertodr/talk-acs)

Browse slides at [http://tinyurl.com/talk-orlando](http://tinyurl.com/talk-orlando)

[remark]: https://github.com/gnab/remark
[cicero]: https://github.com/bast/cicero
