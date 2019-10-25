using FileIO
using Makie

export start_gui, it_exists

picture = 0 # load("./pictures/flowers.jpg")
marker_size = 20
scene_width = 200
scene_height = 200

origin = Point2f0(0, 0)
xybounds = Point2f0(scene_width, scene_height)

function start_gui()
    clicks = Node(Point2f0[])
    scene = scatter(clicks, color = :black, marker = picture,
            markersize = marker_size, limits = FRect(0, 0, xybounds),
            resolution = (500,500))


    on(scene.events.mousebuttons) do buttons
       if ispressed(scene, Mouse.left)
           pos = to_world(scene, Point2f0(scene.events.mouseposition[]))
           if (is_inside_window(pos) & !(it_exists(clicks, pos)))
               push!(clicks, push!(clicks[], pos))
           else
               #println(pos, clicks)
           end
       end
       if ispressed(scene, Mouse.right)
           pos = to_world(scene, Point2f0[scene.events.mouseposition[]])
       end
       return
    end

    # check if the component on that location exists
    function it_exists(clicks::Node, pos::Point2f0)
        valbool = false
        for val in clicks.val
            valbool |= (all(val .< (pos .+ marker_size/2)) .&
                    all(val .> (pos .- marker_size/2)))
        end
        valbool
    end

    # check if pressed location is inside the window
    function is_inside_window(pos)
        return !any(((pos .+ marker_size/2) .> xybounds)
                    .| ((pos .- marker_size/2) .< origin))
    end


    RecordEvents(scene, "mouse_picking")
end
