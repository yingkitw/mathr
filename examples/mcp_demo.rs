//! MCP server demonstration.
//!
//! Shows how to drive mathr's engine the way an AI agent would: through the
//! Model Context Protocol over stdio (JSON-RPC 2.0).
//!
//! Run with:
//!
//! ```text
//! cargo run --example mcp_demo
//! ```
//!
//! Then, from another terminal:
//!
//! ```text
//! cargo run -- mcp
//! ```
//!
//! and paste the JSON-RPC messages this example prints.

use mathr::mcp::{handle_message, parse_json};

fn main() {
    // A realistic agent session: handshake, discover tools, then call a few.
    let messages = [
        // 1. Initialize — the client announces its protocol version.
        r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"mcp-demo","version":"1.0"}}}"#,
        // 2. Acknowledge — a notification, so the server does not reply.
        r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#,
        // 3. Discover the tool catalogue.
        r#"{"jsonrpc":"2.0","id":2,"method":"tools/list"}"#,
        // 4. Evaluate exactly.
        r#"{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"evaluate","arguments":{"expression":"1/2 + 3/4"}}}"#,
        // 5. Differentiate with an explicit variable.
        r#"{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"differentiate","arguments":{"expression":"y*x^2","variable":"y"}}}"#,
        // 6. Integrate, pulling out the free factor `a`.
        r#"{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"integrate","arguments":{"expression":"a*x","variable":"x"}}}"#,
        // 7. A failing tool is a *successful* response with isError — not a
        //    protocol error.
        r#"{"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"evaluate","arguments":{"expression":"@@@"}}}"#,
    ];

    println!("== Messages to send to `mathr mcp` ==\n");
    for m in &messages {
        println!("{m}");
    }

    println!("\n== Responses ==\n");
    for m in &messages {
        let parsed = parse_json(m).expect("valid JSON");
        match handle_message(&parsed) {
            // `None` means the message was a notification — no reply.
            None => println!("(no response — notification)"),
            Some(Ok(response)) => println!("{response}"),
            Some(Err(err)) => println!("error: {}", err),
        }
    }

    // The tool catalogue, one name per line.
    println!("\n== Available tools ==\n");
    for tool in mathr::mcp::tools() {
        println!("  {:<14} {}", tool.name, tool.description);
    }
}
