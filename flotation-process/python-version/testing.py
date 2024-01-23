import keyboard

def on_key_event(e):
    if e.event_type == keyboard.KEY_DOWN:
        print(f"Key {e.name} was pressed")

keyboard.hook(on_key_event)
keyboard.wait("esc")