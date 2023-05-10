tag=`git rev-parse --short HEAD`
container_name="predengine" # change this 
echo BUILDING: sinaiiidgst/${container_name}:$tag
docker build -t sinaiiidgst/${container_name}:${tag} .
docker push sinaiiidgst/${container_name}:${tag}