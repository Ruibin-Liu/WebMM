self.onmessage = e => { if (e.data.type === "ping") self.postMessage({ type: "pong" }); };
