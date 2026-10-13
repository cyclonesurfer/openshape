# ADR 0021 — An agent interface open to any AI agent

Status: accepted, 2026-10-13 (roadmap: first slice in Phase 5, full
command set with the automation API in Phase 9).

## Context

The project wants the application to be drivable by whatever AI agent a
user chooses: an agent should be able to read a document, create and edit
features and sketches, measure, look at the model and export, the way a
user does. The product must not depend on one AI vendor, must stay
GPL-3.0-or-later and Linux-only, and must not let an agent do anything a
user could not see or undo.

The architecture already helps: every change to a document is an
operation in the op log (`ModelOp`), replayed deterministically and
undoable; the UI is one client of that command layer.

## Options

* **A. Embed an AI vendor's SDK** and call its models from inside the
  app. Ties the product to that vendor (and its keys, prices and terms);
  each further vendor is another integration. Rejected.
* **B. Our own JSON-RPC protocol.** Neutral, but no agent speaks it, so
  every agent needs an adapter written for it.
* **C. The Model Context Protocol (MCP)**: an open, published
  specification for exposing tools and resources to agents, implemented by
  most agent clients from different vendors and by open-source agents, and
  governed since December 2025 by the Linux Foundation's Agentic AI
  Foundation. The app is an MCP *server*; the user's agent is the client.
  The app contains no model and no vendor code, and holds no API keys.
* **D. A command-line interface** (`opendatum` subcommands reading and
  writing documents, printing JSON). Any agent that can run a shell can
  use it, MCP or not; it is also what scripts and CI use.

## Decision

**C and D, over one command layer.** The app exposes the same commands
two ways: an MCP server (on standard input and output for a headless
session, and on a local socket to drive the running app, where the user
watches the changes happen), and a CLI with JSON output for agents and
scripts without MCP. Both are thin translations of one typed command
layer, which the automation API of Phase 9 will also use.

Rules:

* **No special powers.** Every agent change is a `ModelOp` through the op
  log: validated like a user's edit, shown in the feature list and
  history, attributed to the agent, undoable. Queries are read-only.
  Nothing an agent does bypasses the kernel's validation.
* **Everything the UI can do is a command.** New UI actions are built as
  commands first (as edits already are, through the op log), so the agent
  interface never lags the UI.
* **Agents need to see.** Besides queries (feature tree, parameters,
  persistent names, measurements, mass properties, errors), a command
  renders a view of the model to an image, using the offscreen renderer
  the tests already use.
* **Neutral and local.** The protocol is implemented from its public
  specification (or with an SDK whose licence `cargo deny` accepts); the
  socket is local to the user's machine and off unless the user turns it
  on; no telemetry.

## Consequences

* The roadmap gains an "Agent interface" step: a first slice in Phase 5
  (documents and history give agents open, save, versions and attributed
  history), with read queries, feature and sketch creation and editing,
  rendering and export; the full command set comes with the automation
  API in Phase 9.
* Agent support costs little while every action stays a command: keeping
  it that way is the one rule that applies before the step itself.
