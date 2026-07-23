#!/usr/bin/env python3

import sys
import re
import html
import base64
from io import BytesIO
import colorsys
import argparse
import gzip
import bz2
import lzma
import datetime
import json

try:
    import compression.zstd
    has_zstd = True
except ModuleNotFoundError:
    has_zstd = False

try:
    from PIL import Image, ImageDraw, ImageFont
except ModuleNotFoundError:
    print(r"""This script requires the Pillow package.
You can install it in a venv with:
    python3 -m venv /tmp/venv/
    /tmp/venv/bin/pip install Pillow
    /tmp/venv/bin/python3 scripts/chronograms/render.py \
            /tmp/chrono.zst -o /tmp/chrono.html
        """)


def feed(files):
    if not files:
        for c in sys.stdin.readlines():
            yield c
    else:
        for file in files:
            if re.match(r"^.*.zstd?$", file):
                assert has_zstd
                with compression.zstd.open(file, 'r') as f:
                    for c in f.readlines():
                        yield c.decode()
            elif re.match(r"^.*.gz$", file):
                with gzip.open(file, 'r') as f:
                    for c in f.readlines():
                        yield c.decode()
            elif re.match(r"^.*.xz$", file):
                with lzma.open(file, 'r') as f:
                    for c in f.readlines():
                        yield c.decode()
            elif re.match(r"^.*.bz2$", file):
                with bz2.open(file) as f:
                    for c in f.readlines():
                        yield c.decode()
            else:
                with open(file) as f:
                    for c in f.readlines():
                        yield c


def convert_epoch_to_readable(t):
    dt = datetime.datetime.fromtimestamp(int(t)).astimezone()
    nsec = int((t - int(t)) * 1e9)
    iso_fmt = f"{dt.strftime('%Y-%m-%d %H:%M:%S')}"
    iso_fmt += f".{nsec:09d}{dt.strftime('%z')}"
    return iso_fmt


# --- Color Palette Management ---
PREDEFINED_COLORS = {
    "INIT": (78, 121, 167),       # Steel Blue
    "SKEWGAUSS": (242, 142, 43),   # Orange
    "ADJUST": (225, 87, 89),      # Soft Red
    "ALLOC": (118, 183, 178),     # Teal
    "FIB": (89, 161, 79),         # Green
    "SSS": (237, 201, 72),        # Gold / Yellow
    "PBR": (176, 122, 161),       # Purple
    "ECM": (255, 157, 167),       # Pink
    "DS": (156, 117, 95),         # Brown
    "PCLAT": (186, 176, 172),     # Warm Gray
    "SLICING": (92, 192, 222),    # Cyan
    "BOTCHED": (200, 0, 0),       # Bright Red
}


def get_color(category):
    """Returns an RGB color tuple for a given event category."""
    cat_upper = category.upper()
    if cat_upper in PREDEFINED_COLORS:
        return PREDEFINED_COLORS[cat_upper]

    hue = (hash(cat_upper) & 0xFFFFFFFF) * 0.618033988749895 % 1.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.85)
    return (int(r * 255), int(g * 255), int(b * 255))


LINE_REGEX = re.compile(r"^t\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+(.+)$")


def parse_trace_data(input_stream):
    bubbles = []
    for line in input_stream:
        line = line.strip()
        if not line or not line.startswith('t '):
            continue
        match = LINE_REGEX.match(line)
        if not match:
            continue

        thread_id = int(match.group(1))
        start_time = float(match.group(2))
        duration = float(match.group(3))
        metric = int(match.group(4))
        description = match.group(5).strip()

        category = description.split()[0] if description else "UNKNOWN"

        bubbles.append({
            "thread": thread_id,
            "start": start_time,
            "duration": duration,
            "end": start_time + duration,
            "metric": metric,
            "desc": description,
            "category": category
        })
    return bubbles


def clamp_trace_data(bubbles, tstart_arg, tend_arg):
    if not bubbles:
        return bubbles, 0.0, 0.0

    raw_min_start = min(b["start"] for b in bubbles)
    raw_max_end = max(b["end"] for b in bubbles)

    if tstart_arg is None:
        win_start = raw_min_start
    elif tstart_arg > 1e8:  # Absolute timestamp
        win_start = tstart_arg
    else:                  # Relative offset
        win_start = raw_min_start + tstart_arg

    if tend_arg is None:
        win_end = raw_max_end
    elif tend_arg > 1e8:    # Absolute timestamp
        win_end = tend_arg
    else:                  # Relative offset
        win_end = raw_min_start + tend_arg

    if win_start >= win_end:
        print(f"Warning: Clamp range [{win_start}, {win_end}] is invalid.",
              file=sys.stderr)
        return [], win_start, win_end

    clamped = []
    for b in bubbles:
        if b["end"] <= win_start or b["start"] >= win_end:
            continue

        cb = dict(b)
        cb["start"] = max(b["start"], win_start)
        cb["end"] = min(b["end"], win_end)
        cb["duration"] = cb["end"] - cb["start"]
        clamped.append(cb)

    return clamped, win_start, win_end


def generate_xhtml_chronogram(bubbles,
                              win_start,
                              win_end,
                              output_filepath="chronogram.html"):
    if not bubbles:
        print("Error: No valid time bubble data found"
              " in the selected time window.", file=sys.stderr)
        return

    threads = sorted(list({b["thread"] for b in bubbles}))
    thread_map = {t: idx for idx, t in enumerate(threads)}
    num_threads = len(threads)

    total_time = max(win_end - win_start, 1e-9)

    # --- Dimensions & Scaling ---
    LEFT_MARGIN = 75
    RIGHT_MARGIN = 20
    HEADER_HEIGHT = 35
    BOTTOM_MARGIN = 25

    target_width = 1200
    per_thread_width = (target_width - LEFT_MARGIN - RIGHT_MARGIN)
    per_thread_width = per_thread_width // num_threads
    per_thread_width = max(4, min(250, per_thread_width))

    img_width = LEFT_MARGIN + RIGHT_MARGIN + num_threads * per_thread_width

    target_height = 12000

    useful_height = target_height - HEADER_HEIGHT - BOTTOM_MARGIN
    default_y_scale = useful_height / total_time

    if False:
        durations = sorted(b["duration"] for b in bubbles)
        p10_idx = int(len(durations) * 0.1)
        target_duration = durations[p10_idx]

        # 2. Compute y_scale so >= 90% of bubbles are at least 2px tall
        if target_duration > 0:
            required_y_scale = 2.0 / target_duration
        else:
            non_zero = [d for d in durations if d > 0]
            required_y_scale = (2.0 / non_zero[0]) if non_zero else default_y_scale

        y_scale = max(default_y_scale, required_y_scale)
        img_height = HEADER_HEIGHT + BOTTOM_MARGIN + int(total_time * y_scale)
        print(f"image height is {img_height}")
    else:
        img_height = target_height
        y_scale = default_y_scale

    img = Image.new("RGB", (img_width, img_height), (255, 255, 255))
    draw = ImageDraw.Draw(img)

    try:
        font = ImageFont.load_default()
    except IOError:
        font = None

    # Draw Grid Lines
    num_ticks = 10
    for i in range(num_ticks + 1):
        rel_t = total_time * (i / num_ticks)
        y = HEADER_HEIGHT + int(rel_t * y_scale)
        draw.line([(LEFT_MARGIN - 5, y),
                   (img_width - RIGHT_MARGIN, y)],
                  fill=(235, 235, 235),
                  width=1)
        draw.text((5, y - 6),
                  f"+{rel_t*1000:.1f}ms",
                  fill=(120, 120, 120),
                  font=font)

    # Draw Column Headers & Dividers
    for t in threads:
        col_idx = thread_map[t]
        x1 = LEFT_MARGIN + col_idx * per_thread_width
        draw.text((x1 + 8, 10), f"Thread {t}", fill=(30, 30, 30), font=font)
        draw.line([(x1, HEADER_HEIGHT),
                   (x1, img_height - BOTTOM_MARGIN)],
                  fill=(210, 210, 210),
                  width=1)

    # Draw Time Bubbles & Build Compact Spatial Index for JS
    trace_data_js = {}

    for b in bubbles:
        col_idx = thread_map[b["thread"]]
        x1 = LEFT_MARGIN + col_idx * per_thread_width
        width_px = per_thread_width

        y1 = HEADER_HEIGHT + int((b["start"] - win_start) * y_scale)
        y_end_calc = HEADER_HEIGHT + int((b["end"] - win_start) * y_scale)

        # Calculate height in pixels
        height_px = max(1, y_end_calc - y1)

        # PIL drawing bounds (inclusive)
        x2_draw = x1 + width_px - 1
        y2_draw = y1 + height_px - 1

        color = get_color(b["category"])
        draw.rectangle([x1, y1, x2_draw, y2_draw], fill=color, outline=None)

        rel_start = b['start'] - win_start
        title_text = (
            f"Thread {b['thread']} | {b['category']} | "
            f"Start: +{rel_start*1000:.3f}ms | "
            f"Duration: {b['duration']*1000:.3f}ms | "
            f"Metric: {b['metric']} | Details: {b['desc']}"
        )

        if col_idx not in trace_data_js:
            trace_data_js[col_idx] = []

        # Store [y1, y2, x1, x2, title_text]
        trace_data_js[col_idx].append([y1, y1 + height_px, x1, x1 + width_px, title_text])

    # Sort each thread's events by y1 to guarantee binary search correctness
    for c_idx in trace_data_js:
        trace_data_js[c_idx].sort(key=lambda x: x[0])

    json_trace_data = json.dumps(trace_data_js, separators=(',', ':'))

    # Encode Base64 Image
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    img_b64 = base64.b64encode(buffer.getvalue()).decode("utf-8")

    categories_present = sorted(list({b["category"] for b in bubbles}))
    legend_html = []
    for cat in categories_present:
        r, g, b = get_color(cat)
        legend_html.append(
            '<span class="legend-item">'
            '<span class="color-box" style="background-color: '
            f'rgb({r},{g},{b});"></span>{cat}</span>'
        )

    xhtml_content = f"""<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html PUBLIC
    "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>cado-nfs chronograms</title>
<style type="text/css">
  html, body {{
    height: 100%;
    margin: 0;
    padding: 0;
    overflow: hidden;
    font-family: -apple-system, BlinkMacSystemFont,
                "Segoe UI", Roboto, Arial, sans-serif;
    background-color: #f0f0f0;
  }}

  .app-container {{
    display: flex;
    flex-direction: column;
    height: 100vh;
    padding: 12px;
    box-sizing: border-box;
  }}

  /* Fixed Header Frame */
  .top-frame {{
    flex: 0 0 auto;
    background: #ffffff;
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 12px 16px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.06);
    margin-bottom: 10px;
  }}

  .header-row {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 6px;
  }}

  .header-title {{
    margin: 0;
    font-size: 16px;
    color: #222;
  }}

  .controls {{
    font-size: 12px;
    color: #444;
    display: flex;
    align-items: center;
    gap: 12px;
    background: #f0f0f0;
    padding: 4px 10px;
    border-radius: 4px;
    border: 1px solid #e0e0e0;
    display: none;
  }}

  .controls label {{
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 3px;
  }}

  .window-meta {{
    font-size: 11px;
    color: #666;
    margin-bottom: 8px;
  }}

  .legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    margin-bottom: 8px;
    padding-bottom: 8px;
    border-bottom: 1px solid #eee;
  }}

  .legend-item {{
    display: inline-flex;
    align-items: center;
    font-size: 12px;
    color: #444;
  }}

  .color-box {{
    width: 12px;
    height: 12px;
    display: inline-block;
    margin-right: 5px;
    border: 1px solid rgba(0,0,0,0.15);
  }}

  #information {{
    margin: 0;
    padding: 8px 12px;
    background: #f8f9fa;
    color: #333333;
    border: 1px solid #dcdcdc;
    border-radius: 3px;
    font-family: SFMono-Regular, Consolas, "Liberation Mono", Menlo, monospace;
    font-size: 12px;
    min-height: 18px;
    white-space: pre-wrap;
    word-break: break-all;
  }}

  /* Scrollable Image Container */
  .chart-scroll-frame {{
    flex: 1 1 auto;
    overflow: auto;
    background: #e9e9e9;
    border: 1px solid #bbb;
    border-radius: 4px;
    position: relative;
  }}

  .image-map-wrapper {{
    position: relative;
    display: inline-block;
    background: #ffffff;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
  }}
</style>
</head>

<body>

<div class="app-container">

  <!-- Fixed Top Panel -->
  <div class="top-frame">
    <div class="header-row">
      <h3 class="header-title">cado-nfs Execution Chronogram</h3>

      <!-- Border Display Radio Option -->
      <div class="controls" id="border-controls">
        <strong>Delimiters:</strong>
        <label><input type="radio" name="border-opt" value="none"
        checked="checked" /> Off</label>
        <label><input type="radio" name="border-opt" value="top" /> On</label>
      </div>
    </div>

    <div class="window-meta">
      Time Window:
      <strong>{convert_epoch_to_readable(win_start)}</strong> to
      <strong>{convert_epoch_to_readable(win_end)}</strong>
      (Span: {(win_end-win_start)*1000:.3f} ms)
      | Rendered Events: {len(bubbles)}
    </div>

    <div class="legend">
      {' '.join(legend_html)}
    </div>

<div id="information">Loading trace data ({len(bubbles)} elements)... </div>
  </div>

  <!-- Scrollable Graphic Container -->
  <div class="chart-scroll-frame">
    <div class="image-map-wrapper">
      <img src="data:image/png;base64,{img_b64}" alt="Chronogram Gantt Chart" />

      <!-- Single Canvas Overlay for Delimiters and Hover Highlights -->
      <canvas id="overlay-canvas"
        width="{img_width}"
        height="{img_height}"
        style="position:absolute; top:0; left:0; cursor:crosshair;">
      </canvas>
    </div>
  </div>

</div>

<script type="text/javascript">
//<![CDATA[
const TRACE = {json_trace_data};
const LEFT_MARGIN = {LEFT_MARGIN};
const PER_THREAD_WIDTH = {per_thread_width};

document.addEventListener("DOMContentLoaded", () => {{
  const infoElement = document.getElementById("information");
  const canvas = document.getElementById("overlay-canvas");
  const ctx = canvas.getContext("2d");
  const controlsElement = document.getElementById("border-controls");

  let showDelimiters = false;
  let activeBubble = null;

  // Binary search within a pre-sorted thread array (O(log N))
  function findBubble(events, mouseY) {{
    if (!events) return null;
    let low = 0, high = events.length - 1;
    while (low <= high) {{
      const mid = (low + high) >> 1;
      const e = events[mid];
      if (mouseY >= e[0] && mouseY < e[1]) {{
        return e;
      }} else if (mouseY < e[0]) {{
        high = mid - 1;
      }} else {{
        low = mid + 1;
      }}
    }}
    return null;
  }}

  function renderOverlay() {{
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Draw top line delimiters if enabled
    if (showDelimiters) {{
      ctx.fillStyle = "rgba(0, 0, 0, 0.6)";
      for (const col in TRACE) {{
        const events = TRACE[col];
        for (let i = 0; i < events.length; i++) {{
          const e = events[i];
          ctx.fillRect(e[2], e[0], e[3] - e[2], 1);
        }}
      }}
    }}

    // Draw hover highlight box
    if (activeBubble) {{
      const [y1, y2, x1, x2] = activeBubble;
      ctx.fillStyle = "rgba(0, 0, 0, 0.35)";
      ctx.fillRect(x1, y1, x2 - x1, y2 - y1);
      ctx.strokeStyle = "#000000";
      ctx.lineWidth = 1;
      ctx.strokeRect(x1, y1, x2 - x1, y2 - y1);
    }}
  }}

  // Radio button border control listener
  document.querySelectorAll("input[name='border-opt']").forEach(radio => {{
    radio.addEventListener("change", (e) => {{
      showDelimiters = (e.target.value === "top");
      renderOverlay();
    }});
  }});

  // Spatial lookup on mouse move (O(1) column index + O(log N) binary search)
  canvas.addEventListener("mousemove", (e) => {{
    const rect = canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    const colIdx = Math.floor((mouseX - LEFT_MARGIN) / PER_THREAD_WIDTH);
    const hit = findBubble(TRACE[colIdx], mouseY);

    if (hit !== activeBubble) {{
      activeBubble = hit;
      if (hit) {{
        infoElement.textContent = hit[4];
      }} else {{
        infoElement.textContent = "Hover over a region to inspect details...";
      }}
      renderOverlay();
    }}
  }});

  canvas.addEventListener("mouseleave", () => {{
    activeBubble = null;
    infoElement.textContent = "Hover over a region to inspect details...";
    renderOverlay();
  }});

  requestAnimationFrame(() => {{
    setTimeout(() => {{
      infoElement.textContent = "Hover over a region to inspect details...";
      controlsElement.style.display = "flex";
    }}, 50);
  }});
}});
//]]>
</script>

</body>
</html>
"""

    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(xhtml_content)

    print(f"Generated Chronogram: {output_filepath}")
    print(f"Clamped Range: {convert_epoch_to_readable(win_start)}"
          f" -> {convert_epoch_to_readable(win_end)}"
          f" (Span: {(win_end-win_start)*1000:.3f} ms)"
          f" ({len(bubbles)} events rendered)")


def main():
    parser = argparse.ArgumentParser(
            description="Generate a Gantt HTML Chronogram"
                        " from CADO-NFS trace logs.")
    parser.add_argument("file", nargs='*',
                        help="Path to trace file (default: stdin)")
    parser.add_argument("-o", "--output", default="chronogram.html",
                        help="Output HTML file path"
                             " (default: chronogram.html)")
    parser.add_argument("--tstart", type=float, default=None,
                        help="Start time offset in seconds"
                             " (relative to run start)"
                             " or absolute timestamp.")
    parser.add_argument("--tend", type=float, default=None,
                        help="End time offset in seconds"
                             " (relative to run start)"
                             " or absolute timestamp.")

    args = parser.parse_args()

    raw_data = parse_trace_data(feed(args.file))
    clamped_data, win_start, win_end = clamp_trace_data(raw_data,
                                                        args.tstart,
                                                        args.tend)

    generate_xhtml_chronogram(clamped_data, win_start, win_end, args.output)


if __name__ == "__main__":
    main()
